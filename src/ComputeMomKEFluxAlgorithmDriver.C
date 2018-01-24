/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include "ComputeMomKEFluxAlgorithmDriver.h"
#include "Algorithm.h"
#include "AlgorithmDriver.h"
#include "FieldTypeDef.h"
#include "Realm.h"
#include "SolutionOptions.h"
#include "master_element/MasterElement.h"

// stk_mesh/base/fem
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

// stk_util
#include <stk_util/parallel/ParallelReduce.hpp>

namespace sierra{
namespace nalu{

class Realm;

//==========================================================================
// Class Definition
//==========================================================================
// ComputeMomKEFluxAlgorithmDriver - Computes Momentum and Kinetic Energy flux balance through the domain
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ComputeMomKEFluxAlgorithmDriver::ComputeMomKEFluxAlgorithmDriver(
  Realm &realm)
  : AlgorithmDriver(realm),
    solnOpts_(*realm.solutionOptions_),
    lumpedMass_(true)
{
    stk::mesh::MetaData & metaData = realm_.meta_data();
    const int nDim = metaData.spatial_dimension();

    solnOpts_.momAlgAccumulation_.resize(nDim);
    solnOpts_.momAlgInflow_.resize(nDim);
    solnOpts_.momAlgOpen_.resize(nDim);
    solnOpts_.momAlgPressureInflow_.resize(nDim);
    solnOpts_.momAlgPressureSymmetry_.resize(nDim); //Has to be zero
    solnOpts_.momAlgPressureWall_.resize(nDim);
    solnOpts_.momAlgPressureOpen_.resize(nDim);
    solnOpts_.momAlgTauInflow_.resize(nDim);
    solnOpts_.momAlgTauSymmetry_.resize(nDim);
    solnOpts_.momAlgTauWall_.resize(nDim);
    solnOpts_.momAlgTauOpen_.resize(nDim);
    solnOpts_.momAlgActSource_.resize(nDim);
    solnOpts_.momAlgTotMom_.resize(nDim);            
    
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
ComputeMomKEFluxAlgorithmDriver::~ComputeMomKEFluxAlgorithmDriver()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- pre_work --------------------------------------------------------
//--------------------------------------------------------------------------
void
ComputeMomKEFluxAlgorithmDriver::pre_work()
{
    // set post processing to zero
    solnOpts_.keAlgAccumulation_ = 0.0;
    solnOpts_.keAlgInflow_ = 0.0;
    solnOpts_.keAlgOpen_ = 0.0;
    solnOpts_.keAlgPressureInflow_ = 0.0;
    solnOpts_.keAlgPressureSymmetry_ = 0.0;
    solnOpts_.keAlgPressureOpen_ = 0.0;  
    solnOpts_.keAlgTauInflow_ = 0.0;
    solnOpts_.keAlgTauSymmetry_ = 0.0;  
    solnOpts_.keAlgTauOpen_ = 0.0;  
    solnOpts_.keAlgDissipation_ = 0.0;
    solnOpts_.keAlgActSourceWork_ = 0.0;
    solnOpts_.keAlgTotTKE_ = 0.0;          

    stk::mesh::MetaData & metaData = realm_.meta_data();
    const int nDim = metaData.spatial_dimension();
    for(int j=0; j<nDim; j++) {
        solnOpts_.momAlgAccumulation_[j] = 0.0;
        solnOpts_.momAlgInflow_[j]= 0.0;
        solnOpts_.momAlgOpen_[j]= 0.0;
        solnOpts_.momAlgPressureInflow_[j] = 0.0;
        solnOpts_.momAlgPressureSymmetry_[j] = 0.0; //Has to be zero
        solnOpts_.momAlgPressureWall_[j] = 0.0; //Has to be zero
        solnOpts_.momAlgPressureOpen_[j] = 0.0;
        solnOpts_.momAlgTauInflow_[j] = 0.0;
        solnOpts_.momAlgTauSymmetry_[j] = 0.0;
        solnOpts_.momAlgTauWall_[j] = 0.0;        
        solnOpts_.momAlgTauOpen_[j] = 0.0;
        solnOpts_.momAlgActSource_[j] = 0.0;
        solnOpts_.momAlgTotMom_[j] = 0.0;                
    }

}

//--------------------------------------------------------------------------
//-------- post_work -------------------------------------------------------
//--------------------------------------------------------------------------
void
ComputeMomKEFluxAlgorithmDriver::post_work()
{

    stk::mesh::MetaData & metaData = realm_.meta_data();
    const int nDim = metaData.spatial_dimension();

    // compute d(rho * u * u)/dt * scv * dt
    double ke_accumulation = 0.0;
    std::vector<double> mom_accumulation(nDim,0.0);
    compute_accumulation(mom_accumulation, ke_accumulation);
    double ke_dissipation = 0.0;
    ke_dissipation = compute_dissipation();
    double act_source_work = 0.0;
    std::vector<double> act_source_force(nDim,0.0);
    compute_act_source_force_work(act_source_force, act_source_work);
    double tot_tke = 0.0;
    std::vector<double> tot_mom(nDim,0.0);
    compute_tot_mom_ke(tot_mom, tot_tke);
    
    // parallel communicate
    std::vector<double> l_sum, g_sum;
    l_sum.resize((1+nDim)*13);
    g_sum.resize((1+nDim)*13);
    l_sum[0] = ke_accumulation;
    l_sum[1] = solnOpts_.keAlgInflow_;
    l_sum[2] = solnOpts_.keAlgOpen_ ;
    l_sum[3] = solnOpts_.keAlgPressureInflow_;
    l_sum[4] = solnOpts_.keAlgPressureSymmetry_; //Has to be zero
    l_sum[5] = solnOpts_.keAlgPressureOpen_;  
    l_sum[6] = solnOpts_.keAlgTauInflow_;
    l_sum[7] = solnOpts_.keAlgTauSymmetry_;  
    l_sum[8] = solnOpts_.keAlgTauWall_;
    l_sum[9] = solnOpts_.keAlgTauOpen_;  
    l_sum[10] = ke_dissipation;
    l_sum[11] = act_source_work;
    l_sum[12] = tot_tke;      
    
    for(int j=0; j<nDim; j++) {
        l_sum[13+j*13]   = mom_accumulation[j];
        l_sum[13+j*13+1] = solnOpts_.momAlgInflow_[j];
        l_sum[13+j*13+2] = solnOpts_.momAlgOpen_[j];
        l_sum[13+j*13+3] = solnOpts_.momAlgPressureInflow_[j];
        l_sum[13+j*13+4] = solnOpts_.momAlgPressureSymmetry_[j]; 
        l_sum[13+j*13+5] = solnOpts_.momAlgPressureOpen_[j];
        l_sum[13+j*13+6] = solnOpts_.momAlgPressureWall_[j];
        l_sum[13+j*13+7] = solnOpts_.momAlgTauInflow_[j];
        l_sum[13+j*13+8] = solnOpts_.momAlgTauSymmetry_[j];
        l_sum[13+j*13+9] = solnOpts_.momAlgTauWall_[j];
        l_sum[13+j*13+10] = solnOpts_.momAlgTauOpen_[j];
        l_sum[13+j*13+11] = act_source_force[j];
        l_sum[13+j*13+12] = tot_mom[j];                        
    }
    stk::ParallelMachine comm = NaluEnv::self().parallel_comm();
    stk::all_reduce_sum(comm, l_sum.data(), g_sum.data(), (1+nDim)*13);
    
    // set parameters for later usage
    solnOpts_.keAlgAccumulation_ = g_sum[0];
    solnOpts_.keAlgInflow_ =  g_sum[1];
    solnOpts_.keAlgOpen_ =  g_sum[2];
    solnOpts_.keAlgPressureInflow_ = g_sum[3];
    solnOpts_.keAlgPressureSymmetry_ = g_sum[4]; //Has to be zero
    solnOpts_.keAlgPressureOpen_ = g_sum[5];  
    solnOpts_.keAlgTauInflow_ = g_sum[6];
    solnOpts_.keAlgTauSymmetry_ = g_sum[7];  
    solnOpts_.keAlgTauWall_ = g_sum[8];
    solnOpts_.keAlgTauOpen_ = g_sum[9];  
    solnOpts_.keAlgDissipation_ = g_sum[10];
    solnOpts_.keAlgActSourceWork_ = g_sum[11];
    solnOpts_.keAlgTotTKE_ = g_sum[12];
    for(int j=0; j<nDim; j++) {
        solnOpts_.momAlgAccumulation_[j] = g_sum[13+j*13];
        solnOpts_.momAlgInflow_[j] = g_sum[13+j*13+1];
        solnOpts_.momAlgOpen_[j] = g_sum[13+j*13+2];
        solnOpts_.momAlgPressureInflow_[j] = g_sum[13+j*13+3];
        solnOpts_.momAlgPressureSymmetry_[j] = g_sum[13+j*13+4]; 
        solnOpts_.momAlgPressureWall_[j] = g_sum[13+j*13+5];
        solnOpts_.momAlgPressureOpen_[j] = g_sum[13+j*13+6];
        solnOpts_.momAlgTauInflow_[j] = g_sum[13+j*13+7];
        solnOpts_.momAlgTauSymmetry_[j] = g_sum[13+j*13+8];
        solnOpts_.momAlgTauWall_[j] = g_sum[13+j*13+9];
        solnOpts_.momAlgTauOpen_[j] = g_sum[13+j*13+10];
        solnOpts_.momAlgActSource_[j] = g_sum[13+j*13+11];
        solnOpts_.momAlgTotMom_[j] = g_sum[13+j*13+12];                
    }

    provide_output();
}

//--------------------------------------------------------------------------
//-------- compute_accumulation --------------------------------------------
//--------------------------------------------------------------------------
void
ComputeMomKEFluxAlgorithmDriver::compute_accumulation(std::vector<double>& mom_accumulation, double & ke_accumulation)
{
  stk::mesh::MetaData & metaData = realm_.meta_data();

  const double dt = realm_.get_time_step();
  const int nDim = metaData.spatial_dimension();

  // initialize accumulation term to zero
  ke_accumulation = 0.0;
  for (int j=0; j<nDim; j++) mom_accumulation[j] = 0.0;

  // extract time parameters
  const double gamma1 = realm_.get_gamma1();
  const double gamma2 = realm_.get_gamma2();
  const double gamma3 = realm_.get_gamma3(); // gamma3 may be zero

  // extract fields
  ScalarFieldType *density = metaData.get_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "density");
  ScalarFieldType &densityNp1 = density->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType &densityN = density->field_of_state(stk::mesh::StateN);
  ScalarFieldType &densityNm1 = (density->number_of_states() == 2) 
    ? density->field_of_state(stk::mesh::StateN) : density->field_of_state(stk::mesh::StateNM1);
  VectorFieldType *velocity = metaData.get_field<VectorFieldType>(
      stk::topology::NODE_RANK, "velocity");
  NaluEnv::self().naluOutput() << " Number of states in velocity =  " << velocity->number_of_states() << std::endl ;
  
  VectorFieldType &velocityNp1 = velocity->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &velocityN = velocity->field_of_state(stk::mesh::StateN);
  VectorFieldType &velocityNm1 = (velocity->number_of_states() == 2) 
      ? velocity->field_of_state(stk::mesh::StateN) : velocity->field_of_state(stk::mesh::StateNM1);
  VectorFieldType *coordinates = metaData.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, solnOpts_.get_coordinates_name());

  //  required space
  std::vector<double> ws_shape_function;
  std::vector<double> ws_rhoNp1;
  std::vector<double> ws_rhoN;
  std::vector<double> ws_rhoNm1;
  std::vector<double> ws_velNp1;
  std::vector<double> ws_velN;
  std::vector<double> ws_velNm1;
  std::vector<double> ws_coordinates;
  std::vector<double> ws_scv_volume;

  // selector (everywhere density lives, locally owned and active) 
  stk::mesh::Selector s_locally_owned = stk::mesh::selectField(*density)    
    & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& elem_buckets =
    realm_.get_buckets( stk::topology::ELEMENT_RANK, s_locally_owned);
  
  for ( stk::mesh::BucketVector::const_iterator ib = elem_buckets.begin() ;
        ib != elem_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    // extract master element
    MasterElement *meSCV = sierra::nalu::MasterElementRepo::get_volume_master_element(b.topology());

    // extract master element specifics
    const int nodesPerElement = meSCV->nodesPerElement_;
    const int numScvIp = meSCV->numIntPoints_;

    // resize
    ws_shape_function.resize(numScvIp*nodesPerElement);
    ws_rhoNp1.resize(nodesPerElement);
    ws_rhoN.resize(nodesPerElement);
    ws_rhoNm1.resize(nodesPerElement);
    ws_velNp1.resize(nodesPerElement*nDim);
    ws_velN.resize(nodesPerElement*nDim);
    ws_velNm1.resize(nodesPerElement*nDim);
    ws_coordinates.resize(nDim*nodesPerElement);
    ws_scv_volume.resize(numScvIp);
    
    if ( lumpedMass_ )
      meSCV->shifted_shape_fcn(&ws_shape_function[0]);
    else
      meSCV->shape_fcn(&ws_shape_function[0]);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      //===============================================
      // gather nodal data; this is how we do it now..
      //===============================================
      stk::mesh::Entity const *  node_rels = b.begin_nodes(k);
      int num_nodes = b.num_nodes(k);

      // sanity check on num nodes
      ThrowAssert( num_nodes == nodesPerElement );

      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = node_rels[ni];

        // pointers to real data
        const double * coords = stk::mesh::field_data(*coordinates, node );

        // gather scalars
        ws_rhoNp1[ni]  = *stk::mesh::field_data(densityNp1, node );
        ws_rhoN[ni]  = *stk::mesh::field_data(densityN, node );
        ws_rhoNm1[ni]  = *stk::mesh::field_data(densityNm1, node );

        // gather vectors
        const int niNdim = ni*nDim;

        const double * uNm1 = stk::mesh::field_data(velocityNm1, node);
        const double * uN   = stk::mesh::field_data(velocityN, node);
        const double * uNp1 = stk::mesh::field_data(velocityNp1, node);
        for ( int j=0; j < nDim; ++j ) {
            ws_velNp1[niNdim+j]  = uNm1[j];
            ws_velN[niNdim+j]  = uN[j];
            ws_velNm1[niNdim+j]  = uNp1[j];
        }
        
        for ( int j=0; j < nDim; ++j ) {
          ws_coordinates[niNdim+j] = coords[j];
        }
      }

      // compute geometry
      double scv_error = 0.0;
      meSCV->determinant(1, &ws_coordinates[0], &ws_scv_volume[0], &scv_error);

      for ( int ip = 0; ip < numScvIp; ++ip ) {

        // zero out; scalar
        double rhoNm1Scv = 0.0;
        double rhoNScv = 0.0;
        double rhoNp1Scv = 0.0;

        std::vector<double> velNm1Scv(nDim,0.0);
        std::vector<double> velNScv(nDim,0.0);
        std::vector<double> velNp1Scv(nDim,0.0);
        
        const int offSet = ip*nodesPerElement;
        for ( int ic = 0; ic < nodesPerElement; ++ic ) {
          // save off shape function
          const double r = ws_shape_function[offSet+ic];
          
          // density
          rhoNm1Scv += r*ws_rhoNm1[ic];
          rhoNScv += r*ws_rhoN[ic];
          rhoNp1Scv += r*ws_rhoNp1[ic];

          for ( int j=0; j < nDim; ++j ) {
              velNm1Scv[j] += r*ws_velNm1[ic*nDim+j];
              velNScv[j] += r*ws_velN[ic*nDim+j];
              velNp1Scv[j] += r*ws_velNp1[ic*nDim+j];              
          }
          
        }

        double keNm1 = 0.0;
        double keN = 0.0;
        double keNp1 = 0.0;        
        
        for (int j=0; j<nDim; j++) {
            keNm1 += 0.5*velNm1Scv[j]*velNm1Scv[j];
            keN += 0.5*velNScv[j]*velNScv[j];
            keNp1 += 0.5*velNp1Scv[j]*velNp1Scv[j];
            mom_accumulation[j] += (gamma1*rhoNp1Scv*velNp1Scv[j] + gamma2*rhoNScv*velNScv[j] + gamma3*rhoNm1Scv*velNm1Scv[j])/dt*ws_scv_volume[ip];
        }

        ke_accumulation +=  
          (gamma1*rhoNp1Scv*keNm1 + gamma2*rhoNScv*keN + gamma3*rhoNm1Scv*keNp1)/dt*ws_scv_volume[ip];
        
      }
    }
  }


  NaluEnv::self().naluOutput() << "  Processor mom accumulation:   " << std::setprecision (6) << mom_accumulation[0] << " " << mom_accumulation[1] << std::endl ;
  NaluEnv::self().naluOutput() << "  Processor ke accumlation:   " << std::setprecision (6) << ke_accumulation << std::endl ;

}

//--------------------------------------------------------------------------
//-------- compute_dissipation ---------------------------------------------
//--------------------------------------------------------------------------
double 
ComputeMomKEFluxAlgorithmDriver::compute_dissipation()
{
  stk::mesh::MetaData & metaData = realm_.meta_data();

  const int nDim = metaData.spatial_dimension();

  // initialize accumulation term to zero
  double ke_dissipation = 0.0;
  
  // extract fields
  GenericFieldType *dudx = metaData.get_field<GenericFieldType>(stk::topology::NODE_RANK, "dudx");

  const std::string viscName = realm_.is_turbulent()
    ? "effective_viscosity_u" : "viscosity";
  ScalarFieldType *viscosity = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, viscName);
  
  VectorFieldType *velocity = metaData.get_field<VectorFieldType>(
      stk::topology::NODE_RANK, "velocity");
  VectorFieldType &velocityNp1 = velocity->field_of_state(stk::mesh::StateNP1);
  VectorFieldType *coordinates = metaData.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, solnOpts_.get_coordinates_name());

  //  required space
  std::vector<double> ws_shape_function;
  std::vector<double> ws_viscosity;
  std::vector<double> ws_velNp1;  
  std::vector<double> ws_dudx;  
  std::vector<double> ws_coordinates;
  std::vector<double> ws_scv_volume;

  // selector (everywhere density lives, locally owned and active) 
  stk::mesh::Selector s_locally_owned = metaData.locally_owned_part() & stk::mesh::selectField(*velocity)    
    & !(realm_.get_inactive_selector() );

  stk::mesh::BucketVector const& elem_buckets =
    realm_.get_buckets( stk::topology::ELEMENT_RANK, s_locally_owned);
  
  for ( stk::mesh::BucketVector::const_iterator ib = elem_buckets.begin() ;
        ib != elem_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    // extract master element
    MasterElement *meSCV = sierra::nalu::MasterElementRepo::get_volume_master_element(b.topology());

    // extract master element specifics
    const int nodesPerElement = meSCV->nodesPerElement_;
    const int numScvIp = meSCV->numIntPoints_;

    // resize
    ws_shape_function.resize(numScvIp*nodesPerElement);
    ws_viscosity.resize(nodesPerElement);
    ws_velNp1.resize(nodesPerElement*nDim);
    ws_dudx.resize(nodesPerElement*nDim*nDim);    
    ws_coordinates.resize(nDim*nodesPerElement);
    ws_scv_volume.resize(numScvIp);
    
    if ( lumpedMass_ )
      meSCV->shifted_shape_fcn(&ws_shape_function[0]);
    else
      meSCV->shape_fcn(&ws_shape_function[0]);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      //===============================================
      // gather nodal data; this is how we do it now..
      //===============================================
      stk::mesh::Entity const *  node_rels = b.begin_nodes(k);
      int num_nodes = b.num_nodes(k);

      // sanity check on num nodes
      ThrowAssert( num_nodes == nodesPerElement );

      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = node_rels[ni];

        // pointers to real data
        const double * coords = stk::mesh::field_data(*coordinates, node );

        // gather vectors
        const int niNdim = ni*nDim;
        const int niNNdim = ni*nDim*nDim;

        const double * uNp1 = stk::mesh::field_data(velocityNp1, node);
        const double * dudxNp1 = stk::mesh::field_data(*dudx, node);
        ws_viscosity[ni] = *stk::mesh::field_data(*viscosity, node);
        for ( int j=0; j < nDim; ++j ) {
            ws_velNp1[niNdim+j]  = uNp1[j];
            for ( int k=0; k < nDim; ++k) {
                ws_dudx[niNNdim+j*nDim+k] = dudxNp1[j*nDim+k];
            }
        }
        
        for ( int j=0; j < nDim; ++j ) {
          ws_coordinates[niNdim+j] = coords[j];
        }
      }

      // compute geometry
      double scv_error = 0.0;
      meSCV->determinant(1, &ws_coordinates[0], &ws_scv_volume[0], &scv_error);

      for ( int ip = 0; ip < numScvIp; ++ip ) {

        // zero out;
        std::vector<double> dudxScv(nDim*nDim,0.0);
        std::vector<double> velNp1Scv(nDim,0.0);
        double viscScv = 0.0;
        
        const int offSet = ip*nodesPerElement;
        for ( int ic = 0; ic < nodesPerElement; ++ic ) {
          // save off shape function
          const double r = ws_shape_function[offSet+ic];
          viscScv += r*ws_viscosity[ic];
          for ( int j=0; j < nDim; ++j ) {
              velNp1Scv[j] += r*ws_velNp1[ic*nDim+j];
              for (int k=0; k < nDim; ++k) {
                  dudxScv[j*nDim+k] += r*ws_dudx[ic*nDim*nDim+j*nDim+k];
              }
          }
        }

        double Pk = 0.0;
        for ( int i = 0; i < nDim; ++i ) {
            const int offSet = nDim*i;
            for ( int j = 0; j < nDim; ++j ) {
                Pk += dudxScv[offSet+j]*(dudxScv[offSet+j] + dudxScv[nDim*j+i]);
            }
        }
        Pk *= viscScv;

        ke_dissipation +=  Pk*ws_scv_volume[ip];
        
      }
    }
  }

  return -ke_dissipation;

}

//--------------------------------------------------------------------------
//-------- compute_act_source_force_work------------------------------------
//--------------------------------------------------------------------------
void 
ComputeMomKEFluxAlgorithmDriver::compute_act_source_force_work(std::vector<double> & act_source_force, double & act_source_work)
{
  stk::mesh::MetaData & metaData = realm_.meta_data();

  const int nDim = metaData.spatial_dimension();

  // initialize accumulation term to zero
  act_source_work = 0.0;
  for (int j=0; j<nDim; j++) act_source_force[j] = 0.0;
  
  // extract fields
  VectorFieldType *actuator_source = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, "actuator_source");

  if (actuator_source != NULL) {

  VectorFieldType *velocity = metaData.get_field<VectorFieldType>(
      stk::topology::NODE_RANK, "velocity");
  VectorFieldType &velocityNp1 = velocity->field_of_state(stk::mesh::StateNP1);
  VectorFieldType *coordinates = metaData.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, solnOpts_.get_coordinates_name());

  //  required space
  std::vector<double> ws_shape_function;
  std::vector<double> ws_viscosity;
  std::vector<double> ws_velNp1;  
  std::vector<double> ws_actuator_source;  
  std::vector<double> ws_coordinates;
  std::vector<double> ws_scv_volume;

  // selector (everywhere density lives, locally owned and active) 
  stk::mesh::Selector s_locally_owned = metaData.locally_owned_part() & stk::mesh::selectField(*velocity)    
    & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& elem_buckets =
    realm_.get_buckets( stk::topology::ELEMENT_RANK, s_locally_owned);
  
  for ( stk::mesh::BucketVector::const_iterator ib = elem_buckets.begin() ;
        ib != elem_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    // extract master element
    MasterElement *meSCV = sierra::nalu::MasterElementRepo::get_volume_master_element(b.topology());

    // extract master element specifics
    const int nodesPerElement = meSCV->nodesPerElement_;
    const int numScvIp = meSCV->numIntPoints_;

    // resize
    ws_shape_function.resize(numScvIp*nodesPerElement);
    ws_viscosity.resize(nodesPerElement);
    ws_velNp1.resize(nodesPerElement*nDim);
    ws_actuator_source.resize(nodesPerElement*nDim);    
    ws_coordinates.resize(nDim*nodesPerElement);
    ws_scv_volume.resize(numScvIp);
    
    if ( lumpedMass_ )
      meSCV->shifted_shape_fcn(&ws_shape_function[0]);
    else
      meSCV->shape_fcn(&ws_shape_function[0]);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      //===============================================
      // gather nodal data; this is how we do it now..
      //===============================================
      stk::mesh::Entity const *  node_rels = b.begin_nodes(k);
      int num_nodes = b.num_nodes(k);

      // sanity check on num nodes
      ThrowAssert( num_nodes == nodesPerElement );

      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = node_rels[ni];

        // pointers to real data
        const double * coords = stk::mesh::field_data(*coordinates, node );

        // gather vectors
        const int niNdim = ni*nDim;

        const double * uNp1 = stk::mesh::field_data(velocityNp1, node);
        const double * actuatorSource= stk::mesh::field_data(*actuator_source, node);
        for ( int j=0; j < nDim; ++j ) {
            ws_velNp1[niNdim+j]  = uNp1[j];
            ws_actuator_source[niNdim+j]  = actuatorSource[j];
            ws_coordinates[niNdim+j] = coords[j];            
        }
      }

      // compute geometry
      double scv_error = 0.0;
      meSCV->determinant(1, &ws_coordinates[0], &ws_scv_volume[0], &scv_error);

      for ( int ip = 0; ip < numScvIp; ++ip ) {

        // zero out;
        std::vector<double> actuatorSourceScv(nDim,0.0);
        std::vector<double> velNp1Scv(nDim,0.0);
        
        const int offSet = ip*nodesPerElement;
        for ( int ic = 0; ic < nodesPerElement; ++ic ) {
          // save off shape function
          const double r = ws_shape_function[offSet+ic];
          for ( int j=0; j < nDim; ++j ) {
              velNp1Scv[j] += r*ws_velNp1[ic*nDim+j];
              actuatorSourceScv[j] += r*ws_actuator_source[ic*nDim+j];
          }
        }

        for ( int i = 0; i < nDim; ++i ) {
            const int offSet = nDim*i;
            act_source_work +=  actuatorSourceScv[i]*ws_scv_volume[ip]*velNp1Scv[i];
            act_source_force[i] += actuatorSourceScv[i]*ws_scv_volume[ip];
        }
      }
    }
  }

  }



}

//--------------------------------------------------------------------------
//-------- compute_tot_mom_ke-----------------------------------------------
//--------------------------------------------------------------------------
void 
ComputeMomKEFluxAlgorithmDriver::compute_tot_mom_ke(std::vector<double> & tot_mom, double & tot_ke)
{
  stk::mesh::MetaData & metaData = realm_.meta_data();

  const int nDim = metaData.spatial_dimension();

  // initialize accumulation term to zero
  tot_ke = 0.0;
  for (int j=0; j<nDim; j++) tot_mom[j] = 0.0;
  
  // extract fields
  ScalarFieldType *density = metaData.get_field<ScalarFieldType>(
      stk::topology::NODE_RANK, "density");
  VectorFieldType *velocity = metaData.get_field<VectorFieldType>(
      stk::topology::NODE_RANK, "velocity");
  VectorFieldType &velocityNp1 = velocity->field_of_state(stk::mesh::StateNP1);
  VectorFieldType *coordinates = metaData.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, solnOpts_.get_coordinates_name());

  //  required space
  std::vector<double> ws_shape_function;
  std::vector<double> ws_density;
  std::vector<double> ws_velNp1;  
  std::vector<double> ws_coordinates;
  std::vector<double> ws_scv_volume;

  // selector (everywhere density lives, locally owned and active) 
  stk::mesh::Selector s_locally_owned = metaData.locally_owned_part() & stk::mesh::selectField(*velocity)    
    & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& elem_buckets =
    realm_.get_buckets( stk::topology::ELEMENT_RANK, s_locally_owned);
  
  for ( stk::mesh::BucketVector::const_iterator ib = elem_buckets.begin() ;
        ib != elem_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    // Extract master element
    MasterElement *meSCV = sierra::nalu::MasterElementRepo::get_volume_master_element(b.topology());

    // extract master element specifics
    const int nodesPerElement = meSCV->nodesPerElement_;
    const int numScvIp = meSCV->numIntPoints_;

    // resize
    ws_shape_function.resize(numScvIp*nodesPerElement);
    ws_density.resize(nodesPerElement);
    ws_velNp1.resize(nodesPerElement*nDim);
    ws_coordinates.resize(nDim*nodesPerElement);
    ws_scv_volume.resize(numScvIp);
    
    if ( lumpedMass_ )
      meSCV->shifted_shape_fcn(&ws_shape_function[0]);
    else
      meSCV->shape_fcn(&ws_shape_function[0]);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      //===============================================
      // gather nodal data; this is how we do it now..
      //===============================================
      stk::mesh::Entity const *  node_rels = b.begin_nodes(k);
      int num_nodes = b.num_nodes(k);

      // sanity check on num nodes
      ThrowAssert( num_nodes == nodesPerElement );

      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = node_rels[ni];

        // pointers to real data
        const double * coords = stk::mesh::field_data(*coordinates, node );

        // gather vectors
        const int niNdim = ni*nDim;

        const double * uNp1 = stk::mesh::field_data(velocityNp1, node);
        ws_density[ni] = *stk::mesh::field_data(*density, node);
        for ( int j=0; j < nDim; ++j ) {
            ws_velNp1[niNdim+j]  = uNp1[j];
            ws_coordinates[niNdim+j] = coords[j];
        }
      }

      // compute geometry
      double scv_error = 0.0;
      meSCV->determinant(1, &ws_coordinates[0], &ws_scv_volume[0], &scv_error);

      for ( int ip = 0; ip < numScvIp; ++ip ) {

        // zero out;
        std::vector<double> velNp1Scv(nDim,0.0);
        double rhoScv = 0.0;
        
        const int offSet = ip*nodesPerElement;
        for ( int ic = 0; ic < nodesPerElement; ++ic ) {
          // save off shape function
          const double r = ws_shape_function[offSet+ic];
          rhoScv += r*ws_density[ic];
          for ( int j=0; j < nDim; ++j ) {
              velNp1Scv[j] += r*ws_velNp1[ic*nDim+j];
          }
        }

        for ( int i = 0; i < nDim; ++i ) {
            const int offSet = nDim*i;
            tot_ke +=  0.5*rhoScv*velNp1Scv[i]*velNp1Scv[i]*ws_scv_volume[ip];
            tot_mom[i] += rhoScv*velNp1Scv[i]*ws_scv_volume[ip];
        }
      }
    }
  }



}


//--------------------------------------------------------------------------
//-------- provide_output -----------------------------------------------
//--------------------------------------------------------------------------
void
ComputeMomKEFluxAlgorithmDriver::provide_output()
{

    stk::mesh::MetaData & metaData = realm_.meta_data();
    const int nDim = metaData.spatial_dimension();
    
    // output momentum closure
    std::vector<double> totalMomClosure(nDim,0.0);
    for(int j=0; j<nDim; j++) {
        totalMomClosure[j] = solnOpts_.momAlgAccumulation_[j] - ( solnOpts_.momAlgInflow_[j] + solnOpts_.momAlgOpen_[j] ) + ( solnOpts_.momAlgPressureInflow_[j] + solnOpts_.momAlgPressureSymmetry_[j] + solnOpts_.momAlgPressureWall_[j] + solnOpts_.momAlgPressureOpen_[j] + solnOpts_.momAlgTauInflow_[j] + solnOpts_.momAlgTauSymmetry_[j] + solnOpts_.momAlgTauWall_[j] + solnOpts_.momAlgTauOpen_[j] + solnOpts_.momAlgActSource_[j] );
    }
    NaluEnv::self().naluOutputP0() << "Momentum Balance Review:  " << std::endl;
    NaluEnv::self().naluOutputP0() << "  Momentum accumulation: " ;
    for(int j=0; j<nDim; j++) NaluEnv::self().naluOutputP0() << solnOpts_.momAlgAccumulation_[j] << " " ;
    NaluEnv::self().naluOutputP0() << std::endl;
    
    NaluEnv::self().naluOutputP0() << "  Total momentum flux:      " ;
    for(int j=0; j<nDim; j++) NaluEnv::self().naluOutputP0() << std::setprecision (16) << solnOpts_.momAlgInflow_[j] + solnOpts_.momAlgOpen_[j]  << " " ;
    NaluEnv::self().naluOutputP0() << std::endl;

    NaluEnv::self().naluOutputP0() << "  Total pressure force:   " ;
    for(int j=0; j<nDim; j++) NaluEnv::self().naluOutputP0() << std::setprecision (6) << solnOpts_.momAlgPressureInflow_[j] + solnOpts_.momAlgPressureSymmetry_[j] + solnOpts_.momAlgPressureWall_[j] + solnOpts_.momAlgPressureOpen_[j] << " " ;
    NaluEnv::self().naluOutputP0() << std::endl;

    NaluEnv::self().naluOutputP0() << "  Total viscous/turbulent force:   " ;
    for(int j=0; j<nDim; j++) NaluEnv::self().naluOutputP0() << std::setprecision (6) << solnOpts_.momAlgTauInflow_[j] + solnOpts_.momAlgTauSymmetry_[j] + solnOpts_.momAlgTauWall_[j] + solnOpts_.momAlgTauOpen_[j] << " " ;
    NaluEnv::self().naluOutputP0() << std::endl;

    NaluEnv::self().naluOutputP0() << "  Total actuator force:   " ;
    for(int j=0; j<nDim; j++) NaluEnv::self().naluOutputP0() << std::setprecision (6) << solnOpts_.momAlgActSource_[j] << " " ;
    NaluEnv::self().naluOutputP0() << std::endl;

    NaluEnv::self().naluOutputP0() << "  Total momentum closure:   " ;
    for(int j=0; j<nDim; j++) NaluEnv::self().naluOutputP0() << std::setprecision (6) << totalMomClosure[j] << " " ;
    NaluEnv::self().naluOutputP0() << std::endl;

    NaluEnv::self().naluOutputP0() << "  Total momentum in domain:   " ;
    for(int j=0; j<nDim; j++) NaluEnv::self().naluOutputP0() << std::setprecision (6) << solnOpts_.momAlgTotMom_[j] << " " ;
    NaluEnv::self().naluOutputP0() << std::endl;
    
    // output kinetic energy closure
    const double totalKEClosure = solnOpts_.keAlgAccumulation_ - ( solnOpts_.keAlgInflow_ + solnOpts_.keAlgOpen_ ) + ( solnOpts_.keAlgPressureInflow_ + solnOpts_.keAlgPressureSymmetry_ + solnOpts_.keAlgPressureOpen_ + solnOpts_.keAlgTauInflow_ + solnOpts_.keAlgTauSymmetry_ + solnOpts_.keAlgTauWall_ + solnOpts_.keAlgTauOpen_ + solnOpts_.keAlgDissipation_ + solnOpts_.keAlgActSourceWork_ );
    NaluEnv::self().naluOutputP0() << "Kinetic Energy Balance Review:  " << std::endl;
    NaluEnv::self().naluOutputP0() << "  Energy accumulation: " << solnOpts_.keAlgAccumulation_ << std::endl;
    NaluEnv::self().naluOutputP0() << "  Total momentum flux:      " << std::setprecision (16) << solnOpts_.keAlgInflow_ + solnOpts_.keAlgOpen_<< std::endl;
    NaluEnv::self().naluOutputP0() << "  Total pressure work:   " << std::setprecision (6) << solnOpts_.keAlgPressureInflow_ + solnOpts_.keAlgPressureSymmetry_ + solnOpts_.keAlgPressureOpen_ << std::endl;
    NaluEnv::self().naluOutputP0() << "  Total viscous work:   " << std::setprecision (6) << solnOpts_.keAlgTauInflow_ + solnOpts_.keAlgTauSymmetry_ + solnOpts_.keAlgTauWall_ + solnOpts_.keAlgTauOpen_ << std::endl;
    NaluEnv::self().naluOutputP0() << "  Viscous/turbulent dissipation:   " << std::setprecision (6) << solnOpts_.keAlgDissipation_ << std::endl;
    NaluEnv::self().naluOutputP0() << "  Actuator Source work:   " << std::setprecision (6) << solnOpts_.keAlgActSourceWork_ << std::endl;
    NaluEnv::self().naluOutputP0() << "  Total kinetic energy closure:   " << std::setprecision (6) << totalKEClosure << std::endl;
    NaluEnv::self().naluOutputP0() << "  Total kinetic energy in domain:   " << std::setprecision (6) << solnOpts_.keAlgTotTKE_ << std::endl;    
    
}

} // namespace nalu
} // namespace Sierra
