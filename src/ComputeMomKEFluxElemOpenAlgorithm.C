/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include "ComputeMomKEFluxElemOpenAlgorithm.h"
#include "Algorithm.h"

#include "FieldTypeDef.h"
#include "Realm.h"
#include "SolutionOptions.h"
#include "master_element/MasterElement.h"

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// ComputeMomKEFluxElemOpenAlgorithm - mdot continuity open bc
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ComputeMomKEFluxElemOpenAlgorithm::ComputeMomKEFluxElemOpenAlgorithm(
  Realm &realm,
  stk::mesh::Part *part)
  : Algorithm(realm, part),
    meshMotion_(realm_.does_mesh_move()),
    velocity_(NULL),
    exposedAreaVec_(NULL),
    openMassFlowRate_(NULL),    
    shiftMomKEFlux_(realm_.get_cvfem_shifted_mdot())
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  velocity_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  exposedAreaVec_ = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "exposed_area_vector");
  openMassFlowRate_ = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "open_mass_flow_rate");
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
ComputeMomKEFluxElemOpenAlgorithm::~ComputeMomKEFluxElemOpenAlgorithm()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
ComputeMomKEFluxElemOpenAlgorithm::execute()
{

  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();
  
  // extract noc
  const std::string dofName = "pressure";
  
  // ip values; both boundary and opposing surface
  std::vector<double> uBip(nDim);

  // pointers to fixed values
  double *p_uBip = &uBip[0];

  // nodal fields to gather
  std::vector<double> ws_vrtm;
  
  // master element
  std::vector<double> ws_shape_function;
  std::vector<double> ws_face_shape_function;

  // deal with interpolation procedure
  const double interpTogether = realm_.get_mdot_interp();
  const double om_interpTogether = 1.0-interpTogether;

  // set accumulation variables
  double keflux_open = 0.0;
  std::vector<double> momflux_open(nDim, 0.0);

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    &stk::mesh::selectUnion(partVec_);

  stk::mesh::BucketVector const& face_buckets =
    realm_.get_buckets( meta_data.side_rank(), s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
        ib != face_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    // face master element
    MasterElement *meFC = sierra::nalu::MasterElementRepo::get_surface_master_element(b.topology());
    const int nodesPerFace = b.topology().num_nodes();
    const int numScsBip = meFC->numIntPoints_;

    // algorithm related; element
    ws_vrtm.resize(nodesPerFace*nDim);
    ws_face_shape_function.resize(numScsBip*nodesPerFace);

    // pointers
    double *p_vrtm = &ws_vrtm[0];
    double *p_face_shape_function = &ws_face_shape_function[0];

    // shape functions; boundary
    if ( shiftMomKEFlux_ )
      meFC->shifted_shape_fcn(&p_face_shape_function[0]);
    else
      meFC->shape_fcn(&p_face_shape_function[0]);
    
    const stk::mesh::Bucket::size_type length   = b.size();

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // get face
      stk::mesh::Entity face = b[k];

      //======================================
      // gather nodal data off of face
      //======================================
      stk::mesh::Entity const * face_node_rels = bulk_data.begin_nodes(face);
      int num_face_nodes = bulk_data.num_nodes(face);
      // sanity check on num nodes
      ThrowAssert( num_face_nodes == nodesPerFace );
      for ( int ni = 0; ni < num_face_nodes; ++ni ) {
        stk::mesh::Entity node = face_node_rels[ni];

        // gather vectors
        double * vrtm = stk::mesh::field_data(*velocity_, node);
        const int offSet = ni*nDim;
        for ( int j=0; j < nDim; ++j ) {
          p_vrtm[offSet+j] = vrtm[j];
        }
      }

      // pointer to face data
      const double * areaVec = stk::mesh::field_data(*exposedAreaVec_, face);
      double * mdot = stk::mesh::field_data(*openMassFlowRate_, face);

      // loop over boundary ips
      for ( int ip = 0; ip < numScsBip; ++ip ) {

        // zero out vector quantities
        for ( int j = 0; j < nDim; ++j ) {
          p_uBip[j] = 0.0;
        }

        // interpolate to bip
        const int offSetSF_face = ip*nodesPerFace;
        for ( int ic = 0; ic < nodesPerFace; ++ic ) {
          const double r = p_face_shape_function[offSetSF_face+ic];
          const int offSetFN = ic*nDim;
          for ( int j = 0; j < nDim; ++j ) {
            p_uBip[j] += r*p_vrtm[offSetFN+j];
          }
        }

        // scatter to mdot and accumulate
        double ke_Bip = 0.0;
        for(int j = 0; j < nDim; ++j) {
            ke_Bip += p_uBip[j] * p_uBip[j] ;
            momflux_open[j] += mdot[ip] * p_uBip[j] ;
        }
        keflux_open += 0.5 * mdot[ip] * ke_Bip;

      }
    }
  }
  // scatter back to solution options; not thread safe
  realm_.solutionOptions_->keAlgOpen_ += keflux_open;
  for(int j = 0; j < nDim; ++j) 
    realm_.solutionOptions_->momAlgOpen_[j] += momflux_open[j];
}

} // namespace nalu
} // namespace Sierra
