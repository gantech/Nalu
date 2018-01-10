/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include "ComputeMomKEFluxEdgeOpenAlgorithm.h"
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
// ComputeMomKEFluxEdgeOpenAlgorithm - compute mdot at edges ip
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ComputeMomKEFluxEdgeOpenAlgorithm::ComputeMomKEFluxEdgeOpenAlgorithm(
  Realm &realm,
  stk::mesh::Part *part)
  : Algorithm(realm, part),
    meshMotion_(realm_.does_mesh_move()),
    velocity_(NULL),
    openMassFlowRate_(NULL)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  velocity_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity_rtm");
  openMassFlowRate_ = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "open_mass_flow_rate");
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
ComputeMomKEFluxEdgeOpenAlgorithm::~ComputeMomKEFluxEdgeOpenAlgorithm()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
ComputeMomKEFluxEdgeOpenAlgorithm::execute()
{
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  // set accumulation variables
  double keflux_open = 0.0;
  std::vector<double> momflux_open(nDim, 0.0);

  // define vector of parent topos; should always be UNITY in size
  std::vector<stk::topology> parentTopo;

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    &stk::mesh::selectUnion(partVec_);

  stk::mesh::BucketVector const& face_buckets =
    realm_.get_buckets( meta_data.side_rank(), s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
        ib != face_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    // extract connected element topology
    b.parent_topology(stk::topology::ELEMENT_RANK, parentTopo);
    ThrowAssert ( parentTopo.size() == 1 );
    stk::topology theElemTopo = parentTopo[0];
    MasterElement *meSCS = sierra::nalu::MasterElementRepo::get_surface_master_element(theElemTopo);

    // size some things that are useful
    const int num_face_nodes = b.topology().num_nodes();
    
    const stk::mesh::Bucket::size_type length   = b.size();

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // pointer to face data
      double * mdot = stk::mesh::field_data(*openMassFlowRate_, b, k);

      // extract the connected element to this exposed face; should be single in size!
      stk::mesh::Entity const * face_elem_rels = b.begin_elements(k);
      ThrowAssert( b.num_elements(k) == 1 );

      // get element; its face ordinal number and populate face_node_ordinals
      stk::mesh::Entity element = face_elem_rels[0];
      const int face_ordinal = b.begin_element_ordinals(k)[0];
      const int *face_node_ordinals = meSCS->side_node_ordinals(face_ordinal);

      // get the relations
      stk::mesh::Entity const * elem_node_rels = bulk_data.begin_nodes(element);

      for ( int ip = 0; ip < num_face_nodes; ++ip ) {

        const int nearestNode = face_node_ordinals[ip];

        // left and right nodes; right is on the face; left is the opposing node
        stk::mesh::Entity nodeR = elem_node_rels[nearestNode];

        const double * vrtm =  stk::mesh::field_data(*velocity_, nodeR );

        // scatter to mdot and accumulate
        double ke_Bip = 0.0;
        for(int j = 0; j < nDim; ++j) {
            ke_Bip += vrtm[j] * vrtm[j] ;
            momflux_open[j] += mdot[ip] * vrtm[j] ;
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
