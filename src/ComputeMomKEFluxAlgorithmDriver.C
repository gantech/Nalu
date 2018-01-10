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

    stk::mesh::MetaData & metaData = realm_.meta_data();
    const int nDim = metaData.spatial_dimension();
    for(int j=0; j<nDim; j++) {
        solnOpts_.momAlgAccumulation_[j] = 0.0;
        solnOpts_.momAlgInflow_[j]= 0.0;
        solnOpts_.momAlgOpen_[j]= 0.0;
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
    double ke_accumulation;
    std::vector<double> mom_accumulation(nDim);
    compute_accumulation(mom_accumulation, ke_accumulation);
    
    // parallel communicate
    std::vector<double> l_sum, g_sum;
    l_sum.resize((1+nDim)*3);
    g_sum.resize((1+nDim)*3);
    l_sum[0] = ke_accumulation;
    l_sum[1] = solnOpts_.keAlgInflow_;
    l_sum[2] = solnOpts_.keAlgOpen_ ;
    for(int j=0; j<nDim; j++) {
        l_sum[3+j*3]   = mom_accumulation[j];
        l_sum[3+j*3+1] = solnOpts_.momAlgInflow_[j];
        l_sum[3+j*3+2] = solnOpts_.momAlgOpen_[j];
    }
    stk::ParallelMachine comm = NaluEnv::self().parallel_comm();
    stk::all_reduce_sum(comm, l_sum.data(), g_sum.data(), (1+nDim)*3);
    
    // set parameters for later usage
    solnOpts_.keAlgAccumulation_ = g_sum[0];
    solnOpts_.keAlgInflow_ =  g_sum[1];
    solnOpts_.keAlgOpen_ =  g_sum[2];
    for(int j=0; j<nDim; j++) {
        solnOpts_.momAlgAccumulation_[j] = g_sum[3+3*nDim];
        solnOpts_.momAlgInflow_[j] = g_sum[3+3*nDim+1];
        solnOpts_.momAlgOpen_[j] = g_sum[3+3*nDim+2];
    }
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
          }
          
        }

        double keNm1 = 0.0;
        double keN = 0.0;
        double keNp1 = 0.0;        
        
        for (int j=0; j<nDim; j++) {
            keNm1 += 0.5*velNm1Scv[j]*velNm1Scv[j];
            keN += 0.5*velNScv[j]*velNScv[j];
            keNp1 += 0.5*velNp1Scv[j]*velNp1Scv[j];
            mom_accumulation[j] = (gamma1*rhoNp1Scv*velNm1Scv[j] + gamma2*rhoNScv*velNm1Scv[j] + gamma3*rhoNm1Scv*velNm1Scv[j])/dt*ws_scv_volume[ip];
        }
        
        ke_accumulation +=  
          (gamma1*rhoNp1Scv*keNm1 + gamma2*rhoNScv*keN + gamma3*rhoNm1Scv*keNp1)/dt*ws_scv_volume[ip];
        
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
        totalMomClosure[j] = solnOpts_.momAlgAccumulation_[j] + solnOpts_.momAlgInflow_[j] + solnOpts_.momAlgOpen_[j];
    }
    NaluEnv::self().naluOutputP0() << "Momentum Balance Review:  " << std::endl;
    NaluEnv::self().naluOutputP0() << "Momentum accumulation: " ;
    for(int j=0; j<nDim; j++) NaluEnv::self().naluOutputP0() << solnOpts_.momAlgAccumulation_[j];
    NaluEnv::self().naluOutputP0() << std::endl;
    
    NaluEnv::self().naluOutputP0() << "Integrated inflow:      " ;
    for(int j=0; j<nDim; j++) NaluEnv::self().naluOutputP0() << std::setprecision (16) << solnOpts_.momAlgInflow_[j];
    NaluEnv::self().naluOutputP0() << std::endl;

    NaluEnv::self().naluOutputP0() << "Integrated open:      " ;
    for(int j=0; j<nDim; j++) NaluEnv::self().naluOutputP0() << std::setprecision (16) << solnOpts_.momAlgOpen_[j];
    NaluEnv::self().naluOutputP0() << std::endl;

    NaluEnv::self().naluOutputP0() << "Total momentum closure:   " ;
    for(int j=0; j<nDim; j++) NaluEnv::self().naluOutputP0() << std::setprecision (6) << totalMomClosure[j] ;
    NaluEnv::self().naluOutputP0() << std::endl;
    
    
  // output kinetic energy closure
    const double totalKEClosure = solnOpts_.keAlgAccumulation_ + solnOpts_.keAlgInflow_ + solnOpts_.keAlgOpen_;
  NaluEnv::self().naluOutputP0() << "Kinetic Energy Balance Review:  " << std::endl;
  NaluEnv::self().naluOutputP0() << "Energy accumulation: " << solnOpts_.keAlgAccumulation_ << std::endl;
  NaluEnv::self().naluOutputP0() << "Integrated inflow:      " << std::setprecision (16) << solnOpts_.keAlgInflow_ << std::endl;
  NaluEnv::self().naluOutputP0() << "Integrated open:      " << std::setprecision (16) << solnOpts_.keAlgOpen_ << std::endl;
  NaluEnv::self().naluOutputP0() << "Total energy closure:   " << std::setprecision (6) << totalKEClosure << std::endl;
  
}

} // namespace nalu
} // namespace Sierra
