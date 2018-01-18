/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include "ComputeMomKEFluxElemWallAlgorithm.h"
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
// ComputeMomKEFluxElemWallAlgorithm - DOCUMENTATION HERE
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ComputeMomKEFluxElemWallAlgorithm::ComputeMomKEFluxElemWallAlgorithm(
  Realm &realm,
  stk::mesh::Part *part)
  : Algorithm(realm, part),
    includeDivU_(realm_.get_divU()),
    meshMotion_(realm_.does_mesh_move()),
    velocity_(NULL),
    exposedAreaVec_(NULL),
    shiftMomKEFlux_(realm_.get_cvfem_shifted_mdot())
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  velocity_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  pressure_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "pressure");
  exposedAreaVec_ = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "exposed_area_vector");
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
ComputeMomKEFluxElemWallAlgorithm::~ComputeMomKEFluxElemWallAlgorithm()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
ComputeMomKEFluxElemWallAlgorithm::execute()
{

  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();
  
  // extract noc
  const std::string dofName = "pressure";
  const bool useShiftedGradOp = realm_.get_shifted_grad_op(dofName);
  std::vector<stk::mesh::Entity> connected_nodes;
  
  // ip values; both boundary and opposing surface
  std::vector<double> uBip(nDim);
  double pressureBip, viscBip ;
  std::vector<double> nx(nDim);
  
  // pointers to fixed values
  double *p_uBip = &uBip[0];
  double *p_nx = &nx[0];

  // nodal fields to gather
  std::vector<double> ws_coordinates;  
  std::vector<double> ws_velocity_face;
  std::vector<double> ws_velocity_elem;  
  std::vector<double> ws_viscosity;
  std::vector<double> ws_pressure;
  
  // master element
  std::vector<double> ws_face_shape_function;
  std::vector<double> ws_dndx;
  std::vector<double> ws_det_j;

  // set accumulation variables
  double keflux_wall = 0.0;
  std::vector<double> momflux_wall(nDim, 0.0);
  double keFlux_PressureWall = 0.0;
  std::vector<double> momFlux_PressureWall(nDim, 0.0);
  double keFlux_TauWall = 0.0;
  std::vector<double> momFlux_TauWall(nDim, 0.0);

  // deal with state
  VectorFieldType &velocityNp1 = velocity_->field_of_state(stk::mesh::StateNP1);

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

    // volume master element
    MasterElement *meSCS = sierra::nalu::MasterElementRepo::get_surface_master_element(theElemTopo);
    const int nodesPerElement = meSCS->nodesPerElement_;
    const int numScsIp = meSCS->numIntPoints_;
    
    // face master element
    MasterElement *meFC = sierra::nalu::MasterElementRepo::get_surface_master_element(b.topology());
    const int nodesPerFace = b.topology().num_nodes();
    const int numScsBip = meFC->numIntPoints_;

    // algorithm related; element
    ws_velocity_elem.resize(nodesPerElement*nDim);
    ws_coordinates.resize(nodesPerElement*nDim);
    ws_dndx.resize(nDim*numScsBip*nodesPerElement);
    ws_det_j.resize(numScsBip);

    ws_velocity_face.resize(nodesPerFace*nDim);
    ws_pressure.resize(nodesPerFace);
    ws_viscosity.resize(nodesPerFace);
    ws_face_shape_function.resize(numScsBip*nodesPerFace);

    // pointers
    double *p_velocity_elem = &ws_velocity_elem[0];    

    double *p_coordinates = &ws_coordinates[0];
    double *p_dndx = &ws_dndx[0];
    double *p_velocity_face = &ws_velocity_face[0];
    double *p_viscosity = &ws_viscosity[0];    
    double *p_pressure = &ws_pressure[0];
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
        double * velocity = stk::mesh::field_data(velocityNp1, node);
        const int offSet = ni*nDim;
        for ( int j=0; j < nDim; ++j ) {
            p_velocity_face[offSet+j] = velocity[j];
        }
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
          p_velocity_elem[offSet+j] = uNp1[j];
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
          
        // zero out vector quantities
        for ( int j = 0; j < nDim; ++j ) {
          p_uBip[j] = 0.0;
        }
        pressureBip = 0.0;
        viscBip = 0.0;
        // interpolate to bip
        for ( int ic = 0; ic < nodesPerFace; ++ic ) {
          const double r = p_face_shape_function[offSetSF_face+ic];
          pressureBip += r*p_pressure[ic];
          viscBip += r*p_viscosity[ic];          
          const int offSetFN = ic*nDim;
          for ( int j = 0; j < nDim; ++j ) {
            p_uBip[j] += r*p_velocity_face[offSetFN+j];
          }
        }

        for(int j = 0; j < nDim; ++j) {
            momFlux_PressureWall[j] += pressureBip  * areaVec[ip*nDim+j] ;
        }

      //================================
      // diffusion second
      //================================
      for ( int ic = 0; ic < nodesPerElement; ++ic ) {
          
          const int offSetDnDx = nDim*nodesPerElement*ip + ic*nDim;
          
          for ( int j = 0; j < nDim; ++j ) {
              
              const double axj = areaVec[faceOffSet+j];
              const double dndxj = p_dndx[offSetDnDx+j];
              const double uxj = p_velocity_elem[ic*nDim+j];
              
              const double divUstress = 2.0/3.0*viscBip*dndxj*uxj*axj*includeDivU_;
              
              for ( int i = 0; i < nDim; ++i ) {
                  
                  // matrix entries
                  int indexR = nearestNode*nDim + i;
                  int rowR = indexR*nodesPerElement*nDim;
                  
                  const double dndxi = p_dndx[offSetDnDx+i];
                  const double uxi = p_velocity_elem[ic*nDim+i];
                  const double nxi = p_nx[i];
                  const double om_nxinxi = 1.0-nxi*nxi;
                 
                  // -mu*dui/dxj*Aj; sneak in divU (explicit)
                  double lhsfac = -viscBip*dndxj*axj;
                  keFlux_TauWall -= (lhsfac*uxi + divUstress)*uxi;
                  momFlux_TauWall[j] -= lhsfac*uxi + divUstress;
                  
                  // -mu*duj/dxi*Aj
                  lhsfac = -viscBip*dndxi*axj;
                  keFlux_TauWall -= (lhsfac*uxj)*uxi;
                  momFlux_TauWall[j] -= lhsfac*uxj ;

              }
          }
      }

      }      
    }
  }
  // scatter back to solution options; not thread safe
  for(int j = 0; j < nDim; ++j) 
      realm_.solutionOptions_->momAlgPressureWall_[j] += momFlux_PressureWall[j];

  realm_.solutionOptions_->keAlgTauWall_ += keFlux_TauWall;
  for(int j = 0; j < nDim; ++j) 
      realm_.solutionOptions_->momAlgTauWall_[j] += momFlux_TauWall[j];
  
  
}
    
} // namespace nalu
} // namespace Sierra
