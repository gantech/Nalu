/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <ComputeMomKEFluxElemSymmetryAlgorithm.h>
#include <EquationSystem.h>
#include <FieldTypeDef.h>
#include <LinearSystem.h>
#include <Realm.h>
#include <SolutionOptions.h>
#include <TimeIntegrator.h>
#include <master_element/MasterElement.h>

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
// ComputeMomKEFluxElemSymmetryAlgorithm - DOCUMENTATION HERE
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ComputeMomKEFluxElemSymmetryAlgorithm::ComputeMomKEFluxElemSymmetryAlgorithm(
  Realm &realm,
  stk::mesh::Part *part)
  : Algorithm(realm, part),
    includeDivU_(realm_.get_divU())
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  velocity_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  pressure_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "pressure");
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  const std::string viscName = realm.is_turbulent()
    ? "effective_viscosity_u" : "viscosity";
  viscosity_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, viscName);
  exposedAreaVec_ = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "exposed_area_vector");
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
ComputeMomKEFluxElemSymmetryAlgorithm::~ComputeMomKEFluxElemSymmetryAlgorithm()
{
    // does nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
ComputeMomKEFluxElemSymmetryAlgorithm::execute()
{

  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();
  
  // extract user options (allow to potentially change over time)
  const std::string dofName = "velocity";
  const bool useShiftedGradOp = realm_.get_shifted_grad_op(dofName);
 
  std::vector<stk::mesh::Entity> connected_nodes;

  // vectors
  std::vector<double> nx(nDim);

  // pointers to fixed values
  double *p_nx = &nx[0];

  // nodal fields to gather
  std::vector<double> ws_velocityNp1;
  std::vector<double> ws_coordinates;
  std::vector<double> ws_viscosity;
  std::vector<double> ws_pressure;
  // master element
  std::vector<double> ws_face_shape_function;
  std::vector<double> ws_dndx;
  std::vector<double> ws_det_j;

  // set accumulation variables
  double keflux_symmetry = 0.0;
  std::vector<double> momflux_symmetry(nDim, 0.0);
  double keFlux_PressureSymmetry = 0.0;
  std::vector<double> momFlux_PressureSymmetry(nDim, 0.0);
  double keFlux_TauSymmetry = 0.0;
  std::vector<double> momFlux_TauSymmetry(nDim, 0.0);

  // deal with state
  VectorFieldType &velocityNp1 = velocity_->field_of_state(stk::mesh::StateNP1);

  // define vector of parent topos; should always be UNITY in size
  std::vector<stk::topology> parentTopo;

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    &stk::mesh::selectUnion(partVec_)
    &!(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& face_buckets =
    realm_.get_buckets( meta_data.side_rank(), s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
        ib != face_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    // extract connected element topology
    b.parent_topology(stk::topology::ELEMENT_RANK, parentTopo);
    ThrowAssert ( parentTopo.size() == 1 );
    stk::topology theElemTopo = parentTopo[0];

    // volume master element
    MasterElement *meSCS = sierra::nalu::MasterElementRepo::get_surface_master_element(theElemTopo);
    const int nodesPerElement = meSCS->nodesPerElement_;

    // face master element
    MasterElement *meFC = sierra::nalu::MasterElementRepo::get_surface_master_element(b.topology());
    const int nodesPerFace = meFC->nodesPerElement_;
    const int numScsBip = meFC->numIntPoints_;

    // resize some things; matrix related
    connected_nodes.resize(nodesPerElement);

    // algorithm related; element
    ws_velocityNp1.resize(nodesPerElement*nDim);
    ws_coordinates.resize(nodesPerElement*nDim);
    ws_viscosity.resize(nodesPerFace);
    ws_pressure.resize(nodesPerFace);
    ws_face_shape_function.resize(numScsBip*nodesPerFace);
    ws_dndx.resize(nDim*numScsBip*nodesPerElement);
    ws_det_j.resize(numScsBip);

    // pointers
    double *p_velocityNp1 = &ws_velocityNp1[0];
    double *p_coordinates = &ws_coordinates[0];
    double *p_viscosity = &ws_viscosity[0];
    double *p_pressure = &ws_pressure[0];
    double *p_face_shape_function = &ws_face_shape_function[0];
    double *p_dndx = &ws_dndx[0];

    // shape function
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
        // gather scalars
        p_pressure[ni] = *stk::mesh::field_data(*pressure_, node);
        p_viscosity[ni] = *stk::mesh::field_data(*viscosity_, node);
      }

      // pointer to face data
      const double * areaVec = stk::mesh::field_data(*exposedAreaVec_, face);

      // extract the connected element to this exposed face; should be single in size!
      stk::mesh::Entity const * face_elem_rels = bulk_data.begin_elements(face);
      ThrowAssert( bulk_data.num_elements(face) == 1 );

      // get element; its face ordinal number
      stk::mesh::Entity element = face_elem_rels[0];
      const stk::mesh::ConnectivityOrdinal* face_elem_ords = bulk_data.begin_element_ordinals(face);
      const int face_ordinal = face_elem_ords[0];

      // mapping from ip to nodes for this ordinal
      const int *ipNodeMap = meSCS->ipNodeMap(face_ordinal);

      //==========================================
      // gather nodal data off of element
      //==========================================
      stk::mesh::Entity const * elem_node_rels = bulk_data.begin_nodes(element);
      int num_nodes = bulk_data.num_nodes(element);
      // sanity check on num nodes
      ThrowAssert( num_nodes == nodesPerElement );
      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = elem_node_rels[ni];
        // set connected nodes
        connected_nodes[ni] = node;
        // gather vectors
        double * uNp1 = stk::mesh::field_data(velocityNp1, node);
        double * coords = stk::mesh::field_data(*coordinates_, node);
        const int offSet = ni*nDim;
        for ( int j=0; j < nDim; ++j ) {
          p_velocityNp1[offSet+j] = uNp1[j];
          p_coordinates[offSet+j] = coords[j];
        }
      }

      // compute dndx
      double scs_error = 0.0;
      if ( useShiftedGradOp )
        meSCS->shifted_face_grad_op(1, face_ordinal, &p_coordinates[0], &p_dndx[0], &ws_det_j[0], &scs_error);
      else
        meSCS->face_grad_op(1, face_ordinal, &p_coordinates[0], &p_dndx[0], &ws_det_j[0], &scs_error);
      
      // loop over boundary ips
      for ( int ip = 0; ip < numScsBip; ++ip ) {

        const int nearestNode = ipNodeMap[ip];

        // offset for bip area vector and types of shape function
        const int faceOffSet = ip*nDim;
        const int offSetSF_face = ip*nodesPerFace;

        // form unit normal
        double asq = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          const double axj = areaVec[faceOffSet+j];
          asq += axj*axj;
        }
        const double amag = std::sqrt(asq);
        for ( int i = 0; i < nDim; ++i ) {
          p_nx[i] = areaVec[faceOffSet+i]/amag;
        }

        // interpolate to bip
        double viscBip = 0.0;
        double pressureBip = 0.0;
        for ( int ic = 0; ic < nodesPerFace; ++ic ) {
          const double r = p_face_shape_function[offSetSF_face+ic];
          viscBip += r*p_viscosity[ic];
          pressureBip += r*p_pressure[ic];          
        }

        for(int j = 0; j < nDim; ++j) {
            momFlux_PressureSymmetry[j] += pressureBip  * areaVec[ip*nDim+j] ;
        }
        
        //================================
        // diffusion second
        //================================
        for ( int ic = 0; ic < nodesPerElement; ++ic ) {

          const int offSetDnDx = nDim*nodesPerElement*ip + ic*nDim;

          for ( int j = 0; j < nDim; ++j ) {

            const double axj = areaVec[faceOffSet+j];
            const double dndxj = p_dndx[offSetDnDx+j];
            const double uxj = p_velocityNp1[ic*nDim+j];

            const double divUstress = 2.0/3.0*viscBip*dndxj*uxj*axj*includeDivU_;

            for ( int i = 0; i < nDim; ++i ) {

              // matrix entries
              int indexR = nearestNode*nDim + i;
              int rowR = indexR*nodesPerElement*nDim;

              const double dndxi = p_dndx[offSetDnDx+i];
              const double uxi = p_velocityNp1[ic*nDim+i];
              const double nxi = p_nx[i];
              const double nxinxi = nxi*nxi;

              // -mu*dui/dxj*Aj*ni*ni; sneak in divU (explicit)
              double lhsfac = -viscBip*dndxj*axj*nxinxi;
              keFlux_TauSymmetry -= (lhsfac*uxi + divUstress*nxinxi)*uxi;
              momFlux_TauSymmetry[j] -= lhsfac*uxi + divUstress*nxinxi;
              
              // -mu*duj/dxi*Aj*ni*ni
              lhsfac = -viscBip*dndxi*axj*nxinxi;
              keFlux_TauSymmetry -= (lhsfac*uxj)*uxi;
              momFlux_TauSymmetry[j] -= lhsfac*uxj ;

            }
          }
        }
      }

      // scatter back to solution options; not thread safe
      realm_.solutionOptions_->keAlgTauSymmetry_ += keFlux_TauSymmetry;
      for(int j = 0; j < nDim; ++j) 
          realm_.solutionOptions_->momAlgTauSymmetry_[j] += momFlux_TauSymmetry[j];

      for(int j = 0; j < nDim; ++j) 
          realm_.solutionOptions_->momAlgPressureSymmetry_[j] += momFlux_PressureSymmetry[j];
      
    }
  }
}

} // namespace nalu
} // namespace Sierra
