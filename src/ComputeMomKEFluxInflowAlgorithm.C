/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include "ComputeMomKEFluxInflowAlgorithm.h"
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
// ComputeMomKEFluxInflowAlgorithm - mdot continuity inflow bc
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ComputeMomKEFluxInflowAlgorithm::ComputeMomKEFluxInflowAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  bool useShifted)
  : Algorithm(realm, part),
    includeDivU_(realm_.get_divU()),
    useShifted_(useShifted),
    velocity_(NULL),
    density_(NULL),
    exposedAreaVec_(NULL)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  velocity_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  const std::string viscName = realm.is_turbulent()
      ? "effective_viscosity_u" : "viscosity";
  viscosity_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, viscName);
  pressure_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "pressure");
  // variable density will need density as a function of user inflow conditions
  density_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  exposedAreaVec_ = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "exposed_area_vector");
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
ComputeMomKEFluxInflowAlgorithm::~ComputeMomKEFluxInflowAlgorithm()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
ComputeMomKEFluxInflowAlgorithm::execute()
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();

  const int nDim = meta_data.spatial_dimension();

  std::vector<stk::mesh::Entity> connected_nodes;
  
  // deal with interpolation procedure
  const double interpTogether = realm_.get_mdot_interp();
  const double om_interpTogether = 1.0-interpTogether;

  // set accumulation variables
  double keFluxInflow = 0.0;
  std::vector<double> momFluxInflow(nDim, 0.0);
  double keFlux_PressureInflow = 0.0;
  std::vector<double> momFlux_PressureInflow(nDim, 0.0);
  double keFlux_TauInflow = 0.0;
  std::vector<double> momFlux_TauInflow(nDim, 0.0);

  // nodal fields to gather; gather everything other than what we are assembling
  std::vector<double> ws_pressure;
  std::vector<double> ws_viscosity;
  std::vector<double> ws_density;
  std::vector<double> ws_velocity_face;
  std::vector<double> ws_velocity_elem;  

  // geometry related to populate
  std::vector<double> ws_coordinates;  
  std::vector<double> ws_shape_function;
  std::vector<double> ws_face_shape_function;
  std::vector<double> ws_dndx;
  std::vector<double> ws_det_j;

  // ip data
  std::vector<double>uBip(nDim);
  std::vector<double>rho_uBip(nDim);
  double pressureBip, viscBip ;
  double *p_uBip = &uBip[0];
  double *p_rho_uBip = &rho_uBip[0];

  ScalarFieldType &densityNp1 = density_->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType &pressureNp1 = pressure_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &velocityNp1 = velocity_->field_of_state(stk::mesh::StateNP1);

  // define vector of parent topos; should always be UNITY in size
  std::vector<stk::topology> parentTopo;
  
  // setup for buckets; union parts and ask for locally owned
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
    
    // extract master element specifics
    MasterElement *meFC = sierra::nalu::MasterElementRepo::get_surface_master_element(b.topology());
    const int nodesPerFace = meFC->nodesPerElement_;
    const int numScsBip = meFC->numIntPoints_;

    // algorithm related
    ws_velocity_elem.resize(nodesPerElement*nDim);
    ws_coordinates.resize(nodesPerElement*nDim);
    connected_nodes.resize(nodesPerElement);
    ws_dndx.resize(nDim*numScsBip*nodesPerElement);
    ws_det_j.resize(numScsBip);
    
    ws_pressure.resize(nodesPerFace);
    ws_viscosity.resize(nodesPerFace);    
    ws_density.resize(nodesPerFace);
    ws_velocity_face.resize(nodesPerFace*nDim);
    ws_shape_function.resize(numScsIp*nodesPerElement);
    ws_face_shape_function.resize(numScsBip*nodesPerFace);

    // pointers
    double *p_velocity_elem = &ws_velocity_elem[0];    
    
    double *p_coordinates = &ws_coordinates[0];
    double *p_dndx = &ws_dndx[0];
    double *p_pressure = &ws_pressure[0];
    double *p_viscosity = &ws_viscosity[0];    
    double *p_density = &ws_density[0];
    double *p_velocity_face = &ws_velocity_face[0];
    double *p_face_shape_function = &ws_face_shape_function[0];

    // shape functions; boundary
    if ( useShifted_ )
      meFC->shifted_shape_fcn(&p_face_shape_function[0]);
    else
      meFC->shape_fcn(&p_face_shape_function[0]);

    const stk::mesh::Bucket::size_type length   = b.size();
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // get face
      stk::mesh::Entity face = b[k];
      
      // face node relations for nodal gather
      stk::mesh::Entity const * face_node_rels = b.begin_nodes(k);

      int num_nodes = b.num_nodes(k);
      for ( int ni = 0; ni < num_nodes; ++ni ) {

        // get the node and form connected_node
        stk::mesh::Entity node = face_node_rels[ni];

        // velocity at nodes
        double * velocity = stk::mesh::field_data(*velocity_, node);
        // gather vectors
        const int offSet = ni*nDim;
        for ( int j=0; j < nDim; ++j ) {
            p_velocity_face[offSet+j] = velocity[j];
        }
        // gather scalar
        p_density[ni] = *stk::mesh::field_data(densityNp1, node);
        p_pressure[ni] = *stk::mesh::field_data(pressureNp1, node);
        p_viscosity[ni] = *stk::mesh::field_data(*viscosity_, node);

      }

      // face data
      double * areaVec = stk::mesh::field_data(*exposedAreaVec_, b, k);

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
      num_nodes = bulk_data.num_nodes(element);
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
      if ( useShifted_ )
        meSCS->shifted_face_grad_op(1, face_ordinal, &p_coordinates[0], &p_dndx[0], &ws_det_j[0], &scs_error);
      else
        meSCS->face_grad_op(1, face_ordinal, &p_coordinates[0], &p_dndx[0], &ws_det_j[0], &scs_error);
      

      for ( int ip = 0; ip < numScsBip; ++ip ) {

        const int nearestNode = ipNodeMap[ip];
          
        // offset for bip area vector and types of shape function
        const int faceOffSet = ip*nDim;
        const int offSetSF_face = ip*nodesPerFace;
        
        // interpolate to scs point; operate on saved off ws_field
        for (int j=0; j < nDim; ++j ) {
          p_uBip[j] = 0.0;
          p_rho_uBip[j] = 0.0;
        }

        double rhoBip = 0.0;
        double pressureBip = 0.0;
        double viscBip = 0.0;
        const int offSet = ip*nodesPerFace;
        for ( int ic = 0; ic < nodesPerFace; ++ic ) {
          const double r = p_face_shape_function[offSet+ic];
          const double rhoIC = p_density[ic];
          const double pressureIC = p_pressure[ic];
          rhoBip += r*p_density[ic];
          viscBip += r*p_viscosity[ic];          
          pressureBip += r*p_pressure[ic];
          for ( int j = 0; j < nDim; ++j ) {
            p_uBip[j] += r*p_velocity_face[ic*nDim+j];
            p_rho_uBip[j] += r*rhoIC*p_velocity_face[ic*nDim+j];
          }
        }

        double mdot = 0.0;
        double ke_Bip = 0.0;
        for ( int j=0; j < nDim; ++j ) {
          mdot += (interpTogether*p_rho_uBip[j] + om_interpTogether*rhoBip*p_uBip[j])*areaVec[ip*nDim+j];
          ke_Bip += p_uBip[j] * p_uBip[j] ;
        }

        for(int j = 0; j < nDim; ++j) {
            momFluxInflow[j] += mdot * p_uBip[j] ;
            momFlux_PressureInflow[j] -= pressureBip  * areaVec[ip*nDim+j] ;
            
        }
        keFluxInflow += 0.5 * mdot * ke_Bip;
        keFlux_PressureInflow -= pressureBip * mdot;

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
                  
                  // -mu*dui/dxj*Aj; sneak in divU (explicit)
                  double lhsfac = -viscBip*dndxj*axj;
                  keFlux_TauInflow -= (lhsfac*uxi + divUstress)*uxi;
                  momFlux_TauInflow[j] -= lhsfac*uxi + divUstress;
                  
                  // -mu*duj/dxi*Aj
                  lhsfac = -viscBip*dndxi*axj;
                  keFlux_TauInflow -= (lhsfac*uxj)*uxi;
                  momFlux_TauInflow[j] -= lhsfac*uxj ;
                  
              }
          }
      }
        
        
      }
    }
  }
  // scatter back to solution options
  realm_.solutionOptions_->keAlgInflow_ += keFluxInflow;
  for(int j=0; j < nDim; ++j) 
      realm_.solutionOptions_->momAlgInflow_[j] += momFluxInflow[j];

  realm_.solutionOptions_->keAlgPressureInflow_ += keFlux_PressureInflow;
  for(int j = 0; j < nDim; ++j) 
      realm_.solutionOptions_->momAlgPressureInflow_[j] += momFlux_PressureInflow[j];
  
  realm_.solutionOptions_->keAlgTauInflow_ += keFlux_TauInflow;
  for(int j = 0; j < nDim; ++j) 
      realm_.solutionOptions_->momAlgTauInflow_[j] += momFlux_TauInflow[j];

}

} // namespace nalu
} // namespace Sierra
