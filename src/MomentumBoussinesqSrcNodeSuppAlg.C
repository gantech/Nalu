/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <MomentumBoussinesqSrcNodeSuppAlg.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <SolutionOptions.h>
#include <SupplementalAlgorithm.h>

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// MomentumBoussinesqSrcNodeSuppAlg - -rho*beta*(T-Tref) g
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
MomentumBoussinesqSrcNodeSuppAlg::MomentumBoussinesqSrcNodeSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    temperature_(NULL),
    dualNodalVolume_(NULL),
    tRef_(298.0),
    rhoRef_(1.0),
    beta_(1.0),
    nDim_(1)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  coordinates_ = meta_data.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, realm_.get_coordinates_name());
  // temperature_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "temperature");
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");

  // extract user parameters from solution options
  tRef_ = realm_.solutionOptions_->referenceTemperature_;
  rhoRef_ = realm_.solutionOptions_->referenceDensity_;
  beta_ = realm_.solutionOptions_->thermalExpansionCoeff_;
  nDim_ = meta_data.spatial_dimension();
  gravity_.resize(nDim_);
  gravity_ = realm_.solutionOptions_->gravity_;
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
MomentumBoussinesqSrcNodeSuppAlg::setup()
{
  // all set up in constructor
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
MomentumBoussinesqSrcNodeSuppAlg::node_execute(
  double */*lhs*/,
  double *rhs,
  stk::mesh::Entity node)
{
  // no lhs contribution; all rhs source term; density should be constant...
  // const double temperature = *stk::mesh::field_data(*temperature_, node );
  const double* coord = stk::mesh::field_data(*coordinates_, node);
  const double zheight = coord[2];

  double temperature = 300.0;
  // Hack in theta variation
  if (zheight <= 650.0) {
    temperature = 300.0;
  }
  else if (zheight <= 750) {
    temperature = 300.0 + (zheight - 650.0) * 0.08; // 8K per 100m
  }
  else {
    temperature = 308.0 + (zheight - 750.0) * 0.003; // 3K per km
  }

  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node );
  const double fac = -rhoRef_*beta_*(temperature - tRef_)*dualVolume;
  for ( int i = 0; i < nDim_; ++i ) {
    rhs[i] += fac*gravity_[i];
  }
}

} // namespace nalu
} // namespace Sierra
