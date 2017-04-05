/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <ActuatorLineFAST.h>
#include <FieldTypeDef.h>
#include <NaluParsing.h>
#include <NaluEnv.h>
#include <Realm.h>
#include <Simulation.h>

// master elements
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

// stk_util
#include <stk_util/parallel/ParallelReduce.hpp>

// stk_search
#include <stk_search/CoarseSearch.hpp>
#include <stk_search/IdentProc.hpp>

// basic c++
#include <vector>
#include <map>
#include <string>
#include <stdexcept>
#include <cmath>

namespace sierra{
namespace nalu{


#ifndef Contiguous2DArrayHack
#define Contiguous2DArrayHack

// Neat hack from http://stackoverflow.com/questions/21943621/how-to-create-a-contiguous-2d-array-in-c to allocate and deallocate contiguous 2D arrays in C++

  /* double **dPtr = create2DArray<double>(10,10); */
  /* dPtr[0][0] = 10;  // for example */
  /* delete2DArray(dPtr);  // free the memory */

template <typename T> T** create2DArray(unsigned nrows, unsigned ncols) {

  T** ptr = new T*[nrows];  // allocate pointers
  T* pool = new T[nrows*ncols];  // allocate pool
  for (unsigned i = 0; i < nrows; ++i, pool += ncols )
    ptr[i] = pool;
  return ptr;
}

template <typename T> void delete2DArray(T** arr) {

  delete [] arr[0];  // remove the pool
  delete [] arr;     // remove the pointers
}

#endif


//==========================================================================
// Class Definition
//==========================================================================
// ActuatorLineFASTInfo - holds all points in the tower specification
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ActuatorLineFASTInfo::ActuatorLineFASTInfo() 
  : procId_(-1),
    procIdGroup_(-1),
    globTurbId_(0),
    numPoints_(0),
    turbineName_("machine_one"),
    coords_(NULL),
    velocity_(NULL),
    force_(NULL),
    nType_(NULL)
{
  boundingSphereVec_.resize(0);
  actuatorLinePointInfoMap_.clear();
}


//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
ActuatorLineFASTInfo::~ActuatorLineFASTInfo()
{
  // nothing to do
  if (MPI_COMM_NULL != turbineComm_) { // Act on this turbine only if I'm a part of the turbine MPI_Group
    delete2DArray(coords_);
    delete2DArray(velocity_);
    delete2DArray(force_);
    delete [] nType_ ;
    delete [] bestX_ ;
  }

  MPI_Group_free(&turbineGroup_);
  MPI_Comm_free(&turbineComm_);
}


//==========================================================================
// Class Definition
//==========================================================================
// ActuatorLineFASTPointInfo - holds individual points search information
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ActuatorLineFASTPointInfo::ActuatorLineFASTPointInfo(double searchRadius):
  searchRadius_(searchRadius)
{
  // nothing to do
}


//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
ActuatorLineFASTPointInfo::~ActuatorLineFASTPointInfo()
{
  // nothing to do
}


//==========================================================================
// Class Definition
//==========================================================================
// ActuatorLineFAST - assemble source term for subgrin turbine; WIP
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ActuatorLineFAST::ActuatorLineFAST(
  Realm &realm,
  const YAML::Node &node)
  : Actuator(realm, node),
    realm_(realm),
    searchMethod_(stk::search::BOOST_RTREE),
    actuatorLineMotion_(false)
{
  // load the data
  load(node);

  /*
    current WIP prototype
    Design concepts:
     1) First and foremost, elements are ghosted to the owning point rank. This 
        probably should be changed since the number of elements might be larger
        than the number of points. Therefore, ghosting points to elements is probably
        easier. This will remove the parallel sum contributions from ghosted elements.
        time will tell..

     2) There can be many specifications with the number of points and omega processed.

     3) in the end, we fill the map of ActuatorLineFASTPointInfo objects and iterate this guy
        to assemble source terms

     4) at present, fake source terms on simple Gaussian weighting

    actuator:
      search_method: stk_octree
      search_target_part: block_1

      specifications:

        - name: machine_zero
          procNo: 0
          epsilon: [ 2.0, 0.0, 0.0 ]
          turbine_pos: [ 0.0, 0.0, 0.0 ]
          restart_filename: "blah"
          FAST_input_filename: "Test01.fst"
          turb_id:  1
          turbine_name: machine_zero
  */
}


//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
ActuatorLineFAST::~ActuatorLineFAST()
{

  FAST.end(); // Call destructors in FAST_cInterface

  MPI_Group_free(&worldMPIGroup_);

  // delete data probes specifications vector
  for ( size_t k = 0; k < actuatorLineInfo_.size(); ++k )
    delete actuatorLineInfo_[k];
}


//--------------------------------------------------------------------------
//-------- compute_elem_force_given_weight ----------------------------------
//--------------------------------------------------------------------------
void
ActuatorLineFAST::compute_elem_force_given_weight(
  const int &nDim,
  const double &g,
  const double *pointForce,
  double *elemForce)
{
  // Multiply the point force by the weight at this element location.
  for ( int j = 0; j < nDim; ++j )
    elemForce[j] = pointForce[j]*g;
}


//--------------------------------------------------------------------------
//-------- isotropic_Gaussian_projection -----------------------------------
//--------------------------------------------------------------------------
double
ActuatorLineFAST::isotropic_Gaussian_projection(
  const int &nDim,
  const double &dis,
  const Coordinates &epsilon)
{
  // Compute the force projection weight at this location using an
  // isotropic Gaussian.
  double g;
  const double pi = acos(-1.0);
  if ( nDim == 2 )
    g = (1.0 / (pow(epsilon.x_,2.0) * pi)) * exp(-pow((dis/epsilon.x_),2.0));
  else
    g = (1.0 / (pow(epsilon.x_,3.0) * pow(pi,1.5))) * exp(-pow((dis/epsilon.x_),2.0));

//std::cout << "g = " << g << std::endl;

  return g;
}


//--------------------------------------------------------------------------
//-------- load ------------------------------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLineFAST::load(
  const YAML::Node & y_node)
{
  // check for any data probes
  const YAML::Node y_actuatorLine = y_node["actuator"];
  if (y_actuatorLine) {
    NaluEnv::self().naluOutputP0() << "ActuatorLineFAST::load" << std::endl;

    // search specifications
    std::string searchMethodName = "na";
    get_if_present(y_actuatorLine, "search_method", searchMethodName, searchMethodName);
    
    // determine search method for this pair
    if ( searchMethodName == "boost_rtree" )
      searchMethod_ = stk::search::BOOST_RTREE;
    else if ( searchMethodName == "stk_octree" )
      searchMethod_ = stk::search::OCTREE;
    else if ( searchMethodName == "stk_kdtree" )
      searchMethod_ = stk::search::KDTREE;
    else
      NaluEnv::self().naluOutputP0() << "ActuatorLineFAST::search method not declared; will use BOOST_RTREE" << std::endl;

    // extract the set of from target names; each spec is homogeneous in this respect
    const YAML::Node searchTargets = y_actuatorLine["search_target_part"];
    if (searchTargets.Type() == YAML::NodeType::Scalar) {
      searchTargetNames_.resize(1);
      searchTargetNames_[0] = searchTargets.as<std::string>() ;
    }
    else {
      searchTargetNames_.resize(searchTargets.size());
      for (size_t i=0; i < searchTargets.size(); ++i) {
        searchTargetNames_[i] = searchTargets[i].as<std::string>() ;
      }
    }

    try {
      FAST.readInputFile(y_actuatorLine);
    }
    catch( const std::runtime_error & ex) {
    }

    // save off number of towers
    const int nTurbinesGlob = FAST.get_nTurbinesGlob() ;

    // each specification can have multiple machines
    for (size_t iTurb = 0; iTurb < nTurbinesGlob; ++iTurb) {
      const YAML::Node cur_turbine = y_actuatorLine["Turbine"+std::to_string(iTurb)];
      ActuatorLineFASTInfo *actuatorLineInfo = new ActuatorLineFASTInfo();
      actuatorLineInfo_.push_back(actuatorLineInfo);
      
      // name
      const YAML::Node theName = cur_turbine["turbine_name"];
      if ( theName )
	actuatorLineInfo->turbineName_ = theName.as<std::string>() ;
      else
	throw std::runtime_error("ActuatorLineFAST: no name provided");
     
      actuatorLineMotion_ = true;
      
      // Force projection function properties
      const YAML::Node epsilon = cur_turbine["epsilon"];
        if ( epsilon )
          actuatorLineInfo->epsilon_ = epsilon.as<Coordinates>() ;
        else
          throw std::runtime_error("ActuatorLineFAST: lacking epsilon vector");
    }
  }
}


//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLineFAST::setup()
{
  // objective: declare the part, register coordinates; must be before populate_mesh()

  double tStart;
  double tEnd;
  double dt ; 
  FAST.setRestart( realm_.restarted_simulation() ) ;
  tStart = realm_.get_current_time() ;
  FAST.setTstart( tStart ) ;
  dt = realm_.get_time_step_from_file();
  FAST.setDt ( dt ) ;
  if ( realm_.get_is_terminate_based_on_time() ) {
    tEnd = realm_.get_total_sim_time() ;
  } 
  else {
    const int ntEnd = realm_.get_max_time_step_count();
    const int ntStart = realm_.get_time_step_count() ;
    tEnd = (ntEnd - ntStart)*dt + tStart ;
  }
  FAST.setTend( tEnd );
  
  MPI_Comm_group(NaluEnv::self().parallel_comm(), &worldMPIGroup_);

}

//--------------------------------------------------------------------------
//-------- allocateTurbinesToProcs------------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLineFAST::allocateTurbinesToProcs()
{
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  stk::mesh::MetaData & metaData = realm_.meta_data();

  // clear some of the search info
  boundingHubSphereVec_.clear();
  boundingProcBoxVec_.clear();

  const int nDim = metaData.spatial_dimension();

  // set all of the candidate elements in the search target names
  populate_candidate_procs();

  const int nTurbinesGlob = FAST.get_nTurbinesGlob() ;
  for (size_t iTurb = 0; iTurb < nTurbinesGlob; ++iTurb) {

    theKey theIdent(NaluEnv::self().parallel_rank(), NaluEnv::self().parallel_rank());

    // define a point that will hold the hub location
    boundingHubBigSphereVec_.clear();
    searchKeyPair_.clear();
    Point hubPointCoords;
    double hubCoords[3] = {};
    FAST.getHubPos(hubCoords, iTurb);
    for ( int j = 0; j < nDim; ++j )  hubPointCoords[j] = hubCoords[j];
    boundingSphere theSphere( Sphere(hubPointCoords, 1.0), theIdent);
    boundingHubSphereVec_.push_back(theSphere);
    boundingSphere theBigSphere( Sphere(hubPointCoords, 95.0), theIdent);
    boundingHubBigSphereVec_.push_back(theBigSphere);
    stk::search::coarse_search( boundingHubBigSphereVec_, boundingProcBoxVec_, searchMethod_, NaluEnv::self().parallel_comm(), searchKeyPair_, false);

    ActuatorLineFASTInfo *actuatorLineInfo = actuatorLineInfo_[iTurb];
    const int nTurbineGroupProcs = searchKeyPair_.size();
    int * turbineGroupProcs = new int[nTurbineGroupProcs];
    int iProc=0;
    std::vector<std::pair<boundingSphere::second_type, boundingElementBox::second_type> >::const_iterator ii;
    for( ii=searchKeyPair_.begin(); ii!=searchKeyPair_.end(); ++ii ) {
      turbineGroupProcs[iProc] = ii->second.proc();
      NaluEnv::self().naluOutput() << "Turbine: " << iTurb << " box_proc: " << ii->second.proc() << std::endl ;
      iProc++;
    }
    MPI_Group_incl(worldMPIGroup_, nTurbineGroupProcs, turbineGroupProcs, &(actuatorLineInfo->turbineGroup_) );
    MPI_Comm_create_group(NaluEnv::self().parallel_comm(), actuatorLineInfo->turbineGroup_, iTurb, &(actuatorLineInfo->turbineComm_) );
  }

  stk::search::coarse_search(boundingHubSphereVec_, boundingProcBoxVec_, searchMethod_, NaluEnv::self().parallel_comm(), searchKeyPair_, false);
  int iTurb=0;
  std::vector<std::pair<boundingSphere::second_type, boundingElementBox::second_type> >::const_iterator ii;
  for( ii=searchKeyPair_.begin(); ii!=searchKeyPair_.end(); ++ii ) {
    const uint64_t theBox = ii->second.id();
    unsigned theRank = NaluEnv::self().parallel_rank();
    const unsigned pt_proc = ii->first.proc();
    const unsigned box_proc = ii->second.proc();

    NaluEnv::self().naluOutput() << "rank: " << theRank << " pt_proc: " << pt_proc << " box_proc: " << box_proc << std::endl ;
    
    FAST.setTurbineProcNo(iTurb, box_proc);
    iTurb++;
  }  

  

}


//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLineFAST::initialize()
{

  allocateTurbinesToProcs();

  FAST.allocateInputData();

  if ( ! FAST.isDryRun() ) {
    FAST.init() ;
  }

  update(); // Update location of actuator points, ghosting etc.
}

//--------------------------------------------------------------------------
//-------- update---- ------------------------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLineFAST::update()
{
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  stk::mesh::MetaData & metaData = realm_.meta_data();
 
  const int nDim = metaData.spatial_dimension();

  // clear some of the search info
  boundingElementBoxVec_.clear();
  searchKeyPair_.clear();

  // set all of the candidate elements in the search target names
  populate_candidate_elements();
  
  // create the ActuatorLineFASTPointInfo
  create_actuator_line_point_info_map();

  // complete filling in the set of elements connected to the centroid
  complete_search();
}


//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLineFAST::execute()
{
  // meta/bulk data and nDim
  stk::mesh::MetaData & metaData = realm_.meta_data();
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  const int nDim = metaData.spatial_dimension();

  // // extract fields
  // VectorFieldType *coordinates 
  //   = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  // VectorFieldType *velocity = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  // VectorFieldType *actuator_source 
  //   = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, "actuator_source");
  // ScalarFieldType *actuator_source_lhs
  //   = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "actuator_source_lhs");
  // ScalarFieldType *g
  //   = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "g");
  // ScalarFieldType *density
  //   = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density"); 
  // // deal with proper viscosity
  // const std::string viscName = realm_.is_turbulent() ? "effective_viscosity" : "viscosity";
  // ScalarFieldType *viscosity
  //   = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, viscName); 

  // // fixed size scratch
  // std::vector<double> ws_pointGasVelocity(nDim);
  // std::vector<double> ws_elemCentroid(nDim);
  // std::vector<double> ws_pointForce(nDim);
  // std::vector<double> ws_elemForce(nDim);
  // double ws_pointGasDensity;
  // double ws_pointGasViscosity;
  // double ws_pointForceLHS;
  
  // // zero out source term; do this manually since there are custom ghosted entities
  // stk::mesh::Selector s_nodes = stk::mesh::selectField(*actuator_source);
  // stk::mesh::BucketVector const& node_buckets =
  //   realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );
  // for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
  //       ib != node_buckets.end() ; ++ib ) {
  //   stk::mesh::Bucket & b = **ib ;
  //   const stk::mesh::Bucket::size_type length   = b.size();
  //   double * actSrc = stk::mesh::field_data(*actuator_source, b);
  //   double * actSrcLhs = stk::mesh::field_data(*actuator_source_lhs, b);
  //   double * gF = stk::mesh::field_data(*g, b);
  //   for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
  //     actSrcLhs[k] = 0.0;
  //     gF[k] = 0.0;
  //     const int offSet = k*nDim;
  //     for ( int j = 0; j < nDim; ++j ) {
  //       actSrc[offSet+j] = 0.0;
  //     }
  //   }
  // }

  // for ( size_t iTurb = 0; iTurb < actuatorLineInfo_.size(); ++iTurb ) {
    
  //   ActuatorLineFASTInfo *ali = actuatorLineInfo_[iTurb];

  //   if (MPI_COMM_NULL != ali->turbineComm) { // Act on this turbine only if I'm a part of the turbine MPI_Group

  //     // loop over map and get velocity at points
  //     std::map<size_t, ActuatorLineFASTPointInfo *>::iterator iterPoint;
  //     for (iterPoint  = ali->actuatorLinePointInfoMap_.begin();
  // 	   iterPoint != ali->actuatorLinePointInfoMap_.end();
  // 	   ++iterPoint) {

  // 	// actuator line info object of interest
  // 	ActuatorLineFASTPointInfo * infoObject = (*iterPoint).second;
  // 	size_t pointId = (*iterPoint).first

  // 	//==========================================================================
  // 	// extract the best element; compute drag given this velocity, property, etc
  // 	// this point drag value will be used by all other elements below
  // 	//==========================================================================
  // 	stk::mesh::Entity bestElem = infoObject->bestElem_;
  // 	if ( bulkData.is_valid(elem) ) {

  // 	  int nodesPerElement = bulkData.num_nodes(bestElem);
	  
  // 	  // resize some work vectors
  // 	  resize_std_vector(nDim, ws_coordinates_, bestElem, bulkData);
  // 	  resize_std_vector(nDim, ws_velocity_, bestElem, bulkData);
  // 	  resize_std_vector(1, ws_viscosity_, bestElem, bulkData);
  // 	  resize_std_vector(1, ws_density_, bestElem, bulkData);
	  
  // 	  // gather nodal data to element nodes; both vector and scalar; coords are used in determinant calc
  // 	  gather_field(nDim, &ws_coordinates_[0], *coordinates, bulkData.begin_nodes(bestElem), 
  // 		       nodesPerElement);
  // 	  gather_field_for_interp(nDim, &ws_velocity_[0], *velocity, bulkData.begin_nodes(bestElem), 
  // 				  nodesPerElement);
  // 	  gather_field_for_interp(1, &ws_viscosity_[0], *viscosity, bulkData.begin_nodes(bestElem), 
  // 				  nodesPerElement);
  // 	  gather_field_for_interp(1, &ws_density_[0], *density, bulkData.begin_nodes(bestElem), 
  // 				  nodesPerElement);
	  
  // 	  // compute volume
  // 	  double elemVolume = compute_volume(nDim, bestElem, bulkData);
	  
  // 	  // interpolate velocity
  // 	  interpolate_field(nDim, bestElem, bulkData, &(infoObject->isoParCoords_[0]), 
  // 			    &ws_velocity_[0], &ws_pointGasVelocity[0]);
	  
  // 	  // interpolate viscosity
  // 	  interpolate_field(1, bestElem, bulkData, &(infoObject->isoParCoords_[0]), 
  // 			    &ws_viscosity_[0], &ws_pointGasViscosity);
	  
  // 	  // interpolate density
  // 	  interpolate_field(1, bestElem, bulkData, &(infoObject->isoParCoords_[0]), 
  // 			    &ws_density_[0], &ws_pointGasDensity);
	  
	  
  // 	  for (int j=0; j<nDim; j++) ali->velocity[pointId][j] = ws_pointGasVelocity_[j];

  // 	} else {
  // 	  for (int j=0; j<nDim; j++) ali->velocity[pointId][j] = 0.0;
  // 	}

  //     }

  //     MPI_Reduce(MPI_IN_PLACE, ali->velocity_[0], ali->numPoints_ * nDim, MPI_DOUBLE, MPI_SUM, ali->procIdGroup_, ali->turbineComm_);

  //     if (ali->procId == NaluEnv::self().parallel_rank()) {
  // 	for (int iNode =0; iNode < ali->numPoints_; iNode++) {
  // 	  if (FAST.isDebug() ) {
  // 	    NaluEnv::self().naluOutput() << "Node " << np << " Velocity = " << ws_pointGasVelocity[0] << " " << ws_pointGasVelocity[1] << " " << ws_pointGasVelocity[2] << " " << std::endl ;
  // 	  }
  // 	  FAST.setVelocity(ws_pointGasVelocity, pointId, infoObject->globTurbId_);
  // 	}
	
      
  // 	if ( ! FAST.isDryRun() ) {
	
  // 	  if ( FAST.isTimeZero() ) {
  // 	    FAST.solution0();
  // 	  }
	  
  // 	  //Step FAST
  // 	  FAST.step();
  // 	}
  //     }

  //   }
  // }
 
  // // Are the actuator points moving?
  // if ( actuatorLineMotion_ )
  //   update();

  // for ( size_t iTurb = 0; iTurb < actuatorLineInfo_.size(); ++iTurb ) {
    
  //   ActuatorLineFASTInfo *ali = actuatorLineInfo_[iTurb];

  //   if (MPI_COMM_NULL != ali->turbineComm) { // Act on this turbine only if I'm a part of the turbine MPI_Group

  //     std::map<size_t, ActuatorLineFASTPointInfo *>::iterator iterPoint;


  //     if ( ali->procId_ == NaluEnv::self().parallel_rank() ) {
  // 	for (iterPoint  = ali->actuatorLinePointInfoMap_.begin();
  // 	     iterPoint != ali->actuatorLinePointInfoMap_.end();
  // 	     ++iterPoint) {

  // 	  // actuator line info object of interest
  // 	  ActuatorLineFASTPointInfo * infoObject = (*iterPoint).second;
  // 	  size_t pointId = (*iterPoint).first ;

  // 	  FAST.getForce(ws_pointForce, pointId, infoObject->globTurbId_);
  // 	  if (FAST.isDebug() ) {
  // 	    NaluEnv::self().naluOutput() << "Node " << pointId << " Type " << ali->nodeType_[pointId] << " Force = " << ws_pointForce[0] << " " << ws_pointForce[1] << " " << ws_pointForce[2] << " " << std::endl ;
  // 	  }
  // 	  for(int j=0; j < nDim; j++) ali->force_[pointId][j] = ws_pointForce[j];
  // 	}
  //     }

  //     for (iterPoint  = ali->actuatorLinePointInfoMap_.begin();
  // 	   iterPoint != ali->actuatorLinePointInfoMap_.end();
  // 	   ++iterPoint) {
	
  // 	// actuator line info object of interest
  // 	ActuatorLineFASTPointInfo * infoObject = (*iterPoint).second;
  // 	size_t pointId = (*iterPoint).first;
  // 	//==========================================================================
  // 	// extract the best element; compute drag given this velocity, property, etc
  // 	// this point drag value will be used by all other elements below
  // 	//==========================================================================
	  
  // 	// get the vector of elements
  // 	std::vector<stk::mesh::Entity> elementVec = infoObject->elementVec_;
	
  // 	// Set up the necessary variables to check that forces/projection function are integrating up correctly.
  // 	double gSum = 0.0;
  // 	double forceSum[nDim];
  // 	forceSum[0] = 0.0;
  // 	forceSum[1] = 0.0;
  // 	if (nDim>2){
  // 	  forceSum[2] = 0.0;
  // 	}
	
  // 	// iterate them and apply source term; gather coords
  // 	for ( size_t k = 0; k < elementVec.size(); ++k ) {
	  
  // 	  stk::mesh::Entity elem = elementVec[k];
	  
  // 	  nodesPerElement = bulkData.num_nodes(elem);
	  
  // 	  // resize some work vectors
  // 	  resize_std_vector(nDim, ws_coordinates_, elem, bulkData);
	  
  // 	  // gather coordinates
  // 	  gather_field(nDim, &ws_coordinates_[0], *coordinates, bulkData.begin_nodes(elem), 
  // 		       nodesPerElement);
	  
  // 	  // compute volume
  // 	  double elemVolume = compute_volume(nDim, elem, bulkData);
	  
  // 	  // determine element centroid
  // 	  compute_elem_centroid(nDim, &ws_elemCentroid[0], nodesPerElement);
	  
  // 	  // compute distance
  // 	  const double distance = compute_distance(nDim, &ws_elemCentroid[0], &(infoObject->centroidCoords_[0]));
	  
  // 	  double gA = 0.0;
	  
  // 	  switch (infoObject->nodeType_) {
  // 	  case HUB:
  // 	    // project the force to this element centroid with projection function
  // 	    gA = isotropic_Gaussian_projection(nDim, distance, infoObject->epsilon_);
  // 	    compute_elem_force_given_weight(nDim, gA, &ws_pointForce[0], &ws_elemForce[0]);
  // 	    // assemble nodal quantity; no LHS contribution here...
  // 	    assemble_source_to_nodes(nDim, elem, bulkData, elemVolume, &ws_elemForce[0], ws_pointForceLHS,
  // 				     0.0, *actuator_source, *actuator_source_lhs, *g, 0.0);
  // 	    forceSum[0] += ws_elemForce[0]*elemVolume;
  // 	    forceSum[1] += ws_elemForce[1]*elemVolume;
  // 	    if (nDim > 2){
  // 	      forceSum[2] += ws_elemForce[2]*elemVolume;
  // 	    }
  // 	    gSum += gA*elemVolume;
  // 	    break;
	    
  // 	  case BLADE:
  // 	    // project the force to this element centroid with projection function
  // 	    gA = isotropic_Gaussian_projection(nDim, distance, infoObject->epsilon_);
  // 	    compute_elem_force_given_weight(nDim, gA, &ws_pointForce[0], &ws_elemForce[0]);
  // 	    // assemble nodal quantity; no LHS contribution here...
  // 	    assemble_source_to_nodes(nDim, elem, bulkData, elemVolume, &ws_elemForce[0], ws_pointForceLHS,
  // 				     gA, *actuator_source, *actuator_source_lhs, *g, 0.0);
  // 	    forceSum[0] += ws_elemForce[0]*elemVolume;
  // 	    forceSum[1] += ws_elemForce[1]*elemVolume;
  // 	    if (nDim > 2){
  // 	      forceSum[2] += ws_elemForce[2]*elemVolume;
  // 	    }
  // 	    gSum += gA*elemVolume;
  // 	    break;
	    
  // 	  case TOWER:
  // 	    // project the force to this element centroid with projection function
  // 	    gA = isotropic_Gaussian_projection(nDim, distance, infoObject->epsilon_);
  // 	    compute_elem_force_given_weight(nDim, gA, &ws_pointForce[0], &ws_elemForce[0]);
  // 	    // assemble nodal quantity; no LHS contribution here...
  // 	    assemble_source_to_nodes(nDim, elem, bulkData, elemVolume, &ws_elemForce[0], ws_pointForceLHS,
  // 				     gA, *actuator_source, *actuator_source_lhs, *g, 0.0);
  // 	    forceSum[0] += ws_elemForce[0]*elemVolume;
  // 	    forceSum[1] += ws_elemForce[1]*elemVolume;
  // 	    if (nDim > 2){
  // 	      forceSum[2] += ws_elemForce[2]*elemVolume;
  // 	    }
  // 	    gSum += gA*elemVolume;
  // 	    break;	
  // 	  }
  // 	}

  // 	NaluEnv::self().naluOutput() << "Actuator Point " << np << ", " << "# elems " << elementVec.size() << ", " << " Type " << infoObject->nodeType_ << std::endl;
  // 	NaluEnv::self().naluOutput() << "  -Body force = " << forceSum[0] << " " << forceSum[1] << " " << forceSum[2] << " " << std::endl;
  // 	NaluEnv::self().naluOutput() << "  -Act. force = " << ws_pointForce[0] << " " << ws_pointForce[1] << " " << ws_pointForce[2] << std::endl;
  // 	NaluEnv::self().naluOutput() << "  -Ratio = " << forceSum[0]/ws_pointForce[0] << " " << forceSum[1]/ws_pointForce[1] << " " << forceSum[2]/ws_pointForce[2] << std::endl;
  // 	NaluEnv::self().naluOutput() << "  -gSum = " << gSum << std::endl;
	
  // 	np=np+1;
  //     }

  //     // parallel assemble (contributions from ghosted and locally owned)
  //     const std::vector<const stk::mesh::FieldBase*> sumFieldVec(1, actuator_source);
  //     stk::mesh::parallel_sum_including_ghosts(bulkData, sumFieldVec);
      
  //     const std::vector<const stk::mesh::FieldBase*> sumFieldG(1, g);
  //     stk::mesh::parallel_sum_including_ghosts(bulkData, sumFieldG);
  //   }
    
  // }

}
//--------------------------------------------------------------------------
//-------- populate_candidate_elements -------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLineFAST::populate_candidate_elements() 
{
  stk::mesh::MetaData & metaData = realm_.meta_data();
  stk::mesh::BulkData & bulkData = realm_.bulk_data();

  const int nDim = metaData.spatial_dimension();

  // fields
  VectorFieldType *coordinates = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  // point data structures
  Point minCorner, maxCorner;

  // extract part
  stk::mesh::PartVector searchParts;
  for ( size_t k = 0; k < searchTargetNames_.size(); ++k ) {
    stk::mesh::Part *thePart = metaData.get_part(searchTargetNames_[k]);
    if ( NULL != thePart )
      searchParts.push_back(thePart);
    else
      throw std::runtime_error("ActuatorLineFAST: Part is null" + searchTargetNames_[k]);     
  }

  // selector and bucket loop
  stk::mesh::Selector s_locally_owned = metaData.locally_owned_part()
    &stk::mesh::selectUnion(searchParts);
  
  stk::mesh::BucketVector const& elem_buckets =
    realm_.get_buckets( stk::topology::ELEMENT_RANK, s_locally_owned );

  for ( stk::mesh::BucketVector::const_iterator ib = elem_buckets.begin();
        ib != elem_buckets.end() ; ++ib ) {
    
    stk::mesh::Bucket & b = **ib;

    const stk::mesh::Bucket::size_type length   = b.size();

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // get element
      stk::mesh::Entity elem = b[k];

      // initialize max and min
      for (int j = 0; j < nDim; ++j ) {
        minCorner[j] = +1.0e16;
        maxCorner[j] = -1.0e16;
      }

      // extract elem_node_relations
      stk::mesh::Entity const* elem_node_rels = bulkData.begin_nodes(elem);
      const int num_nodes = bulkData.num_nodes(elem);

      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = elem_node_rels[ni];
        
        // pointers to real data
        const double * coords = stk::mesh::field_data(*coordinates, node );
        
        // check max/min
        for ( int j = 0; j < nDim; ++j ) {
          minCorner[j] = std::min(minCorner[j], coords[j]);
          maxCorner[j] = std::max(maxCorner[j], coords[j]);
        }
      }
      
      // setup ident
      stk::search::IdentProc<uint64_t,int> theIdent(bulkData.identifier(elem), NaluEnv::self().parallel_rank());
      
      // create the bounding point box and push back
      boundingElementBox theBox(Box(minCorner,maxCorner), theIdent);
      boundingElementBoxVec_.push_back(theBox);
    }
  }
}

//--------------------------------------------------------------------------
//-------- populate_candidate_procs ----------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLineFAST::populate_candidate_procs() 
{
  stk::mesh::MetaData & metaData = realm_.meta_data();
  stk::mesh::BulkData & bulkData = realm_.bulk_data();

  const int nDim = metaData.spatial_dimension();

  // fields
  VectorFieldType *coordinates = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  // point data structures
  //  Point minCorner, maxCorner;
  std::vector<Point> minCorner(1), maxCorner(1);
  std::vector<Point> gMinCorner, gMaxCorner;

  // initialize max and min
  for (int j = 0; j < nDim; ++j ) {
    minCorner[0][j] = +1.0e16;
    maxCorner[0][j] = -1.0e16;
  }

  // extract part
  stk::mesh::PartVector searchParts;
  for ( size_t k = 0; k < searchTargetNames_.size(); ++k ) {
    stk::mesh::Part *thePart = metaData.get_part(searchTargetNames_[k]);
    if ( NULL != thePart )
      searchParts.push_back(thePart);
    else
      throw std::runtime_error("ActuatorLineFAST: Part is null" + searchTargetNames_[k]);     
  }

  // selector and bucket loop
  stk::mesh::Selector s_locally_owned = metaData.locally_owned_part()
    &stk::mesh::selectUnion(searchParts);
  
  stk::mesh::BucketVector const& elem_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_locally_owned );

  for ( stk::mesh::BucketVector::const_iterator ib = elem_buckets.begin();
        ib != elem_buckets.end() ; ++ib ) {
    
    stk::mesh::Bucket & b = **ib;

    const stk::mesh::Bucket::size_type length   = b.size();

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // get element
      stk::mesh::Entity node = b[k];

      // pointers to real data
      const double * coords = stk::mesh::field_data(*coordinates, node );
        
      // check max/min
      for ( int j = 0; j < nDim; ++j ) {
	minCorner[0][j] = std::min(minCorner[0][j], coords[j]);
	maxCorner[0][j] = std::max(maxCorner[0][j], coords[j]);
      }

    }      
  }

  stk::parallel_vector_concat(NaluEnv::self().parallel_comm(), minCorner, gMinCorner);
  stk::parallel_vector_concat(NaluEnv::self().parallel_comm(), maxCorner, gMaxCorner);
  
  for(unsigned j = 0; j < NaluEnv::self().parallel_size(); j++) {
    // setup ident
    stk::search::IdentProc<uint64_t,int> theIdent(j, j);
    NaluEnv::self().naluOutput() << "proc " << j << " minCorner: " << gMinCorner[j] << " maxCorner: " << gMaxCorner[j] << std::endl ;
    
    // create the bounding point box and push back
    boundingElementBox theBox(Box(gMinCorner[j],gMaxCorner[j]), theIdent);
    boundingProcBoxVec_.push_back(theBox);
    
  }

}

//--------------------------------------------------------------------------
//-------- create_actuator_line_point_info_map -----------------------------
//--------------------------------------------------------------------------
void
ActuatorLineFAST::create_actuator_line_point_info_map() {

  const double currentTime = realm_.get_current_time();

  stk::mesh::MetaData & metaData = realm_.meta_data(); 
  const int nDim = metaData.spatial_dimension();

  size_t np = 0;

  for ( size_t iTurb = 0; iTurb < actuatorLineInfo_.size(); ++iTurb ) {
    
    ActuatorLineFASTInfo *ali = actuatorLineInfo_[iTurb];

    if (MPI_COMM_NULL != ali->turbineComm_) { // Act on this turbine only if I'm a part of the turbine MPI_Group

      ali->procId_ = FAST.get_procNo(iTurb); 
      int turbGroupProcId = 0; // id of the processor owning FAST in the turbine group comm
      if ( ali->procId_ == NaluEnv::self().parallel_rank() ) {
  	ali->numPoints_ = FAST.get_numForcePts(iTurb); // Total number of actuator points for the turbine
  	MPI_Comm_rank(ali->turbineComm_, &turbGroupProcId);
      }
      MPI_Allreduce(&turbGroupProcId, &(ali->procIdGroup_), 1, MPI_INT, MPI_SUM, ali->turbineComm_) ; // Set the group processor id owning the turbine on all processors in the group
      MPI_Bcast(&(ali->numPoints_), 1, MPI_INT, ali->procIdGroup_, ali->turbineComm_) ; 

      ali->coords_ = create2DArray<double> (ali->numPoints_,3);
      ali->velocity_ = create2DArray<double> (ali->numPoints_,3);
      ali->force_ = create2DArray<double> (ali->numPoints_,3);
      ali->nType_ = new int[ali->numPoints_];
      ali->bestX_ = new minDist[ali->numPoints_];
      ali->boundingSphereVec_.resize(ali->numPoints_);

      if ( ali->procId_ == NaluEnv::self().parallel_rank() ) {
	
  	if (! FAST.isDryRun() ) {
  	  // scratch array for coordinates
  	  double currentCoords[3] = {};
  	  for(int iNode = 0; iNode < ali->numPoints_; iNode++) {
  	    FAST.getForceNodeCoordinates(currentCoords, np, iTurb);
  	    for (int j=0; j < nDim; j++) ali->coords_[iNode][j] = currentCoords[j] ;
  	    ali->nType_[iNode] = FAST.getVelNodeType(np, iTurb);
  	    if (FAST.isDebug() ) {
  	      NaluEnv::self().naluOutput() << "Vel Node " << np << " Position = " << currentCoords[0] << " " << currentCoords[1] << " " << currentCoords[2] << " " << std::endl ;
  	    }
  	  }
  	}
      }
      
      MPI_Bcast(ali->coords_[0], ali->numPoints_*nDim, MPI_DOUBLE, ali->procIdGroup_, ali->turbineComm_) ; 
      MPI_Bcast(ali->nType_, ali->numPoints_*nDim, MPI_INT, ali->procIdGroup_, ali->turbineComm_) ; 

      for(int iNode = 0; iNode < ali->numPoints_; iNode++) {
  	// define a point that will hold the actuator point
  	Point centroidCoords;
  	for(int j=0; j < nDim; j++) centroidCoords[j] = ali->coords_[iNode][j];

  	stk::search::IdentProc<uint64_t,int> theIdent(iNode, NaluEnv::self().parallel_rank());
  	double searchRadius = ali->epsilon_.x_ * sqrt(log(1.0/0.001));

  	// create the bounding point sphere and push back
  	boundingSphere theSphere( Sphere(centroidCoords, searchRadius), theIdent);
  	ali->boundingSphereVec_.push_back(theSphere);
	
  	// create the point info and push back to map
  	ActuatorLineFASTPointInfo *actuatorLinePointInfo 
  	  = new ActuatorLineFASTPointInfo(searchRadius);
  	ali->actuatorLinePointInfoMap_[iNode] = actuatorLinePointInfo;
  	}

    }
  }
}


//--------------------------------------------------------------------------
//-------- complete_search -------------------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLineFAST::complete_search()
{
  stk::mesh::MetaData & metaData = realm_.meta_data();
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  const int nDim = metaData.spatial_dimension();

  // extract fields
  VectorFieldType *coordinates 
    = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  // for ( size_t iTurb = 0; iTurb < actuatorLineInfo_.size(); ++iTurb ) {
    
  //   ActuatorLineFASTInfo *ali = actuatorLineInfo_[iTurb];

  //   if (MPI_COMM_NULL != ali->turbineComm_) { // Act on this turbine only if I'm a part of the turbine MPI_Group

  //     stk::search::coarse_search(ali->boundingSphereVec_, boundingElementBoxVec_, searchMethod_, 
  // 				 NaluEnv::self().parallel_comm(), searchKeyPair_, false); // Search only within this processor

  //     // now proceed with the standard search
  //     std::vector<std::pair<boundingSphere::second_type, boundingElementBox::second_type> >::const_iterator ii;
  //     for( ii=searchKeyPair_.begin(); ii!=searchKeyPair_.end(); ++ii ) {

  // 	const uint64_t thePt = ii->first.id();
  // 	const uint64_t theBox = ii->second.id();
  // 	const unsigned theRank = NaluEnv::self().parallel_rank();
  // 	const unsigned pt_proc = ii->first.proc();
	

  // 	// proceed as required; all elements should have already been ghosted via the coarse search
  // 	stk::mesh::Entity elem = bulkData.get_entity(stk::topology::ELEMENT_RANK, theBox);
  // 	if ( !(bulkData.is_valid(elem)) )
  // 	  throw std::runtime_error("no valid entry for element");

  // 	// find the point data structure
  // 	std::map<size_t, ActuatorLineFASTPointInfo *>::iterator iterPoint;
  // 	iterPoint=ali->actuatorLinePointInfoMap_.find(thePt);
  // 	if ( iterPoint == actuatorLinePointInfoMap_.end() )
  // 	  throw std::runtime_error("no valid entry for actuatorLinePointInfoMap_");
      
  // 	// extract the point object and push back the element to either the best 
  // 	// candidate or the standard vector of elements
  // 	ActuatorLineFASTPointInfo *actuatorLinePointInfo = iterPoint->second;
	
  // 	// extract topo and master element for this topo
  // 	const stk::mesh::Bucket &theBucket = bulkData.bucket(elem);
  // 	const stk::topology &elemTopo = theBucket.topology();
  // 	MasterElement *meSCS = realm_.get_surface_master_element(elemTopo);
  // 	const int nodesPerElement = meSCS->nodesPerElement_;

  // 	// gather elemental coords
  // 	std::vector<double> elementCoords(nDim*nodesPerElement);
  // 	gather_field(nDim, &elementCoords[0], *coordinates, bulkData.begin_nodes(elem), 
  // 		     nodesPerElement);

  // 	// find isoparametric points
  // 	std::vector<double> isoParCoords(nDim);
  // 	const double nearestDistance = meSCS->isInElement(&elementCoords[0],
  // 							  &(actuatorLinePointInfo->centroidCoords_[0]),
  // 							  &(isoParCoords[0]));
	
  // 	// save off best element and its isoparametric coordinates for this point
  // 	if ( nearestDistance < actuatorLinePointInfo->bestX_ ) {
  // 	  actuatorLinePointInfo->bestX_ = nearestDistance;
  // 	  actuatorLinePointInfo->isoParCoords_ = isoParCoords;
  // 	  actuatorLinePointInfo->bestElem_ = elem;
  // 	}
  //       actuatorLinePointInfo->elementVec_.push_back(elem);
  //     }
  //   }
  // }
  
}

//--------------------------------------------------------------------------
//-------- resize_std_vector -----------------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLineFAST::resize_std_vector( 
  const int &sizeOfField,
  std::vector<double> &theVector,   
  stk::mesh::Entity elem, 
  const stk::mesh::BulkData & bulkData)
{
  const stk::topology &elemTopo = bulkData.bucket(elem).topology();
  MasterElement *meSCS = realm_.get_surface_master_element(elemTopo);
  const int nodesPerElement = meSCS->nodesPerElement_;
  theVector.resize(nodesPerElement*sizeOfField);
}


//--------------------------------------------------------------------------
//-------- gather_field ----------------------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLineFAST::gather_field(
  const int &sizeOfField,
  double *fieldToFill, 
  const stk::mesh::FieldBase &stkField,
  stk::mesh::Entity const* elem_node_rels, 
  const int &nodesPerElement) 
{
  for ( int ni = 0; ni < nodesPerElement; ++ni ) { 
    stk::mesh::Entity node = elem_node_rels[ni];     
    const double * theField = (double*)stk::mesh::field_data(stkField, node );
    for ( int j = 0; j < sizeOfField; ++j ) { 
      const int offSet = ni*sizeOfField+j;
      fieldToFill[offSet] = theField[j];
    }   
  }   
}


//--------------------------------------------------------------------------
//-------- gather_field_for_interp -----------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLineFAST::gather_field_for_interp(
  const int &sizeOfField,
  double *fieldToFill, 
  const stk::mesh::FieldBase &stkField,
  stk::mesh::Entity const* elem_node_rels, 
  const int &nodesPerElement) 
{
  for ( int ni = 0; ni < nodesPerElement; ++ni ) { 
    stk::mesh::Entity node = elem_node_rels[ni];     
    const double * theField = (double*)stk::mesh::field_data(stkField, node );
    for ( int j = 0; j < sizeOfField; ++j ) { 
      const int offSet = j*nodesPerElement + ni; 
      fieldToFill[offSet] = theField[j];
    }   
  }   
}


//--------------------------------------------------------------------------
//-------- compute_volume --------------------------------------------------
//--------------------------------------------------------------------------
double
ActuatorLineFAST::compute_volume(
  const int &nDim,
  stk::mesh::Entity elem, 
  const stk::mesh::BulkData & bulkData) 
{
  // extract master element from the bucket in which the element resides
  const stk::topology &elemTopo = bulkData.bucket(elem).topology();
  MasterElement *meSCV = realm_.get_volume_master_element(elemTopo);
  int nodesPerElement = meSCV->nodesPerElement_;
  const int numScvIp = meSCV->numIntPoints_;

  // compute scv for this element
  ws_scv_volume_.resize(numScvIp);
  double scv_error = 0.0;
  meSCV->determinant(1, &ws_coordinates_[0], &ws_scv_volume_[0], &scv_error);

  double elemVolume = 0.0;
  for ( int ip = 0; ip < numScvIp; ++ip ) {
    elemVolume += ws_scv_volume_[ip];
  }
  return elemVolume;
}


//--------------------------------------------------------------------------
//-------- interpolate_field -----------------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLineFAST::interpolate_field(
  const int &sizeOfField,
  stk::mesh::Entity elem, 
  const stk::mesh::BulkData & bulkData,
  double *isoParCoords,
  const double *fieldAtNodes,
  double *pointField) 
{
  // extract master element from the bucket in which the element resides
  const stk::topology &elemTopo = bulkData.bucket(elem).topology();
  MasterElement *meSCS = realm_.get_surface_master_element(elemTopo);
  
  // interpolate velocity to this best point
  meSCS->interpolatePoint(
    sizeOfField,
    isoParCoords,
    fieldAtNodes,
    pointField); 
}


//--------------------------------------------------------------------------
//-------- compute_elem_centroid -------------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLineFAST::compute_elem_centroid(
  const int &nDim,
  double *elemCentroid,
  const int & nodesPerElement) 
{
  // zero
  for ( int j = 0; j < nDim; ++j )
    elemCentroid[j] = 0.0;
  
  // assemble
  for ( int ni = 0; ni < nodesPerElement; ++ni ) {
    for ( int j=0; j < nDim; ++j ) {
      elemCentroid[j] += ws_coordinates_[ni*nDim+j]/nodesPerElement;
    }
  }
}


//--------------------------------------------------------------------------
//-------- compute_distance ------------------------------------------------
//--------------------------------------------------------------------------
double
ActuatorLineFAST::compute_distance(
  const int &nDim,
  const double *elemCentroid,
  const double *pointCentroid) 
{ 
  double distance = 0.0;
  for ( int j = 0; j < nDim; ++j )
    distance += std::pow(elemCentroid[j] - pointCentroid[j], 2);
  distance = std::sqrt(distance);
  return distance;
}


//--------------------------------------------------------------------------
//-------- assemble_source_to_nodes ----------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLineFAST::assemble_source_to_nodes(
  const int &nDim,
  stk::mesh::Entity elem,
  const stk::mesh::BulkData & bulkData,
  const double &elemVolume,
  const double *drag,
  const double &dragLHS,
  const double &gLocal,
  stk::mesh::FieldBase &actuator_source,
  stk::mesh::FieldBase &actuator_source_lhs,
  stk::mesh::FieldBase &g,
  const double &lhsFac)
{
  // extract master element from the bucket in which the element resides
  const stk::topology &elemTopo = bulkData.bucket(elem).topology();
  MasterElement *meSCV = realm_.get_volume_master_element(elemTopo);
  int nodesPerElement = meSCV->nodesPerElement_;
  const int numScvIp = meSCV->numIntPoints_;

  // extract elem_node_relations
  stk::mesh::Entity const* elem_node_rels = bulkData.begin_nodes(elem);

  // assemble to nodes
  const int *ipNodeMap = meSCV->ipNodeMap();
  for ( int ip = 0; ip < numScvIp; ++ip ) {

    // nearest node to ip
    const int nearestNode = ipNodeMap[ip];

    // extract node and pointer to source term
    stk::mesh::Entity node = elem_node_rels[nearestNode];
    double * sourceTerm = (double*)stk::mesh::field_data(actuator_source, node );
    double * sourceTermLHS = (double*)stk::mesh::field_data(actuator_source_lhs, node );
    double * gGlobal = (double*)stk::mesh::field_data(g, node);


    // nodal weight based on volume weight
    const double nodalWeight = ws_scv_volume_[ip]/elemVolume;
    *sourceTermLHS += nodalWeight*dragLHS*lhsFac;
    *gGlobal += gLocal;
    for ( int j=0; j < nDim; ++j ) {
      sourceTerm[j] += nodalWeight*drag[j];
    }
  }
} 
 

} // namespace nalu
} // namespace Sierra
