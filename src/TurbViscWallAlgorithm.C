/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <TurbViscWallAlgorithm.h>
#include <Algorithm.h>
#include <CopyFieldAlgorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <master_element/MasterElement.h>
#include <ABLProfileFunction.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetEntities.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// TurbViscWallAlgorithm - compute tvisc for Wall model
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
TurbViscWallAlgorithm::TurbViscWallAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  const bool &useShifted,
  const double &gravity,
  const double &z0,
  const double &Tref)
  : Algorithm(realm, part),
    useShifted_(useShifted),
    z0_(z0), 
    Tref_(Tref), 
    gravity_(gravity), 
    alpha_h_(1.0), 
    beta_m_(16.0), 
    beta_h_(16.0),
    gamma_m_(5.0),
    gamma_h_(5.0),
    kappa_(realm.get_turb_model_constant(TM_kappa)),    
    bcHeatFlux_(NULL),
    specificHeat_(NULL),
    density_(NULL),
    velocity_(NULL),
    tvisc_(NULL),
    tviscWall_(NULL),
    assembledWallArea_(NULL),
    wallFrictionVelocityBip_(NULL),
    wallNormalDistanceBip_(NULL),
    dudx_(NULL),
    assembledWallNormalDistance_(NULL)
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  density_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  velocity_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  specificHeat_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "specific_heat");
  bcHeatFlux_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "heat_flux_bc");  
  tvisc_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_viscosity");
  tviscWall_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "tvisc_wall");
  exposedAreaVec_ = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "exposed_area_vector");
  assembledWallArea_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "assembled_wall_area_wf");
  wallFrictionVelocityBip_ = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "wall_friction_velocity_bip");
  wallNormalDistanceBip_ = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "wall_normal_distance_bip");  
  assembledWallNormalDistance_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "assembled_wall_normal_distance");
  dudx_ = meta_data.get_field<GenericFieldType>(stk::topology::NODE_RANK, "dudx");

}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
TurbViscWallAlgorithm::execute()
{

  ABLProfileFunction *p_ABLProfFun;
  StableABLProfileFunction StableProfFun(gamma_m_, gamma_h_);
  UnstableABLProfileFunction UnstableProfFun(beta_m_, beta_h_);
  NeutralABLProfileFunction NeutralProfFun;

  stk::mesh::MetaData & meta_data = realm_.meta_data();
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();

  const int nDim = meta_data.spatial_dimension();
  const double invNdim = 1.0/meta_data.spatial_dimension();

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
      & stk::mesh::selectUnion(partVec_);

  stk::mesh::BucketVector const& node_buckets =
      realm_.get_buckets( stk::topology::NODE_RANK, s_locally_owned_union );

  //Zero nodal field on the wall
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
        ib != node_buckets.end() ; ++ib ) {
      stk::mesh::Bucket & b = **ib ;
      const stk::mesh::Bucket::size_type length   = b.size();
      double *tviscWall = stk::mesh::field_data(*tviscWall_, *b.begin() );
      for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
          tviscWall[k] = 0.0;
      }
  }

  // nodal fields to gather
  std::vector<double> ws_velocityNp1;
  std::vector<double> ws_dudx;
  std::vector<double> ws_bcHeatFlux;
  std::vector<double> ws_density;
  std::vector<double> ws_specificHeat;

  // master element
  std::vector<double> ws_face_shape_function;

  // deal with state
  VectorFieldType &velocityNp1 = velocity_->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType &densityNp1 = density_->field_of_state(stk::mesh::StateNP1);

  // bip values
  std::vector<double> uBip(nDim);
  std::vector<double> dudxBip(nDim*nDim);
  std::vector<double> unitNormal(nDim);

  // pointers to fixed values
  double *p_uBip = &uBip[0];
  double *p_dudxBip = &dudxBip[0];
  double *p_unitNormal= &unitNormal[0];
  
  stk::mesh::BucketVector const& face_buckets =
    realm_.get_buckets( meta_data.side_rank(), s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
        ib != face_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    // face master element
    MasterElement *meFC = sierra::nalu::MasterElementRepo::get_surface_master_element(b.topology());
    const int nodesPerFace = meFC->nodesPerElement_;
    const int numScsBip = meFC->numIntPoints_;

    // mapping from ip to nodes for this ordinal; face perspective (use with face_node_relations)
    const int *faceIpNodeMap = meFC->ipNodeMap();

    // algorithm related; element
    ws_velocityNp1.resize(nodesPerFace*nDim);
    ws_dudx.resize(nodesPerFace*nDim*nDim);
    ws_bcHeatFlux.resize(nodesPerFace);
    ws_density.resize(nodesPerFace);
    ws_specificHeat.resize(nodesPerFace);

    ws_face_shape_function.resize(numScsBip*nodesPerFace);

    // pointers
    double *p_velocityNp1 = &ws_velocityNp1[0];
    double *p_dudx = &ws_dudx[0];
    double *p_bcHeatFlux = &ws_bcHeatFlux[0];
    double *p_density = &ws_density[0];
    double *p_specificHeat = &ws_specificHeat[0];
    double *p_face_shape_function = &ws_face_shape_function[0];

    // shape functions
    if ( useShifted_ )
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
      for ( int ni = 0; ni < nodesPerFace; ++ni ) {
        stk::mesh::Entity node = face_node_rels[ni];

        // gather scalars
        p_bcHeatFlux[ni] = *stk::mesh::field_data(*bcHeatFlux_, node);
        p_density[ni]    = *stk::mesh::field_data(densityNp1, node);
        p_specificHeat[ni] = *stk::mesh::field_data(*specificHeat_, node);
        
        // gather vectors
        double * uNp1 = stk::mesh::field_data(velocityNp1, node);
        const double * dudx = (double*)stk::mesh::field_data(*dudx_, node);
        const int offSet = ni*nDim;
        const int offSetT = ni * nDim * nDim ;
        for ( int j=0; j < nDim; ++j ) {
          p_velocityNp1[offSet+j] = uNp1[j];
          for ( int i=0; i < nDim; ++i) {
              p_dudx[offSetT + nDim*j + i] = dudx[nDim*j + i];
          }
        }
      }

      // pointer to face data
      const double * areaVec = stk::mesh::field_data(*exposedAreaVec_, face);
      const double *wallNormalDistanceBip = stk::mesh::field_data(*wallNormalDistanceBip_, face);
      const double *wallFrictionVelocityBip = stk::mesh::field_data(*wallFrictionVelocityBip_, face);

      // loop over face nodes
      for ( int ip = 0; ip < numScsBip; ++ip ) {

        const int offSetAveraVec = ip*nDim;
        const int offSetSF_face = ip*nodesPerFace;

        const int localFaceNode = faceIpNodeMap[ip];
        stk::mesh::Entity lfNode = face_node_rels[localFaceNode];

        // zero out vector quantities; squeeze in aMag
        double aMag = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          p_uBip[j] = 0.0;
          for ( int i = 0; i < nDim; ++i )
              p_dudxBip[j*nDim+i]= 0.0;
          const double axj = areaVec[offSetAveraVec+j];
          aMag += axj*axj;
        }
        aMag = std::sqrt(aMag);

        // interpolate to bip
        double heatFluxBip = 0.0;
        double rhoBip = 0.0;
        double CpBip = 0.0;
        for ( int ic = 0; ic < nodesPerFace; ++ic ) {
          const double r = p_face_shape_function[offSetSF_face+ic];
          rhoBip += r*p_density[ic];
          CpBip += r*p_specificHeat[ic];
          heatFluxBip += r*p_bcHeatFlux[ic];
          const int offSetFN = ic*nDim;
          const int offSetFNT = ic*nDim*nDim;
          for ( int j = 0; j < nDim; ++j ) {
            p_uBip[j] += r*p_velocityNp1[offSetFN+j];
            for ( int i = 0; i < nDim; ++i ) {
                p_dudxBip[j*nDim + i] += r*p_dudx[offSetFNT+j*nDim+i];
            }
          }
        }

        // form unit normal
        for ( int j = 0; j < nDim; ++j ) {
          p_unitNormal[j] = areaVec[offSetAveraVec+j]/aMag;
        }

        // extract bip data
        const double yp = wallNormalDistanceBip[ip];
        const double utau= wallFrictionVelocityBip[ip];

        // determine lambda = tau_w / (rho * u_* * u_tan)
        // NOTE:
        // density should be obtained from the wall function for
        // temperature and a relation rho = rho(T)
        const double TfluxBip = heatFluxBip / (rhoBip * CpBip);

        const double eps_heat_flux = 1.0e-8;
        const double largenum = 1.0e8;
        double Lfac;
        if (TfluxBip < eps_heat_flux) {
          p_ABLProfFun = &StableProfFun;
        }
        else if (TfluxBip > eps_heat_flux) {
          p_ABLProfFun = &UnstableProfFun;
        }
        else {
          p_ABLProfFun = &NeutralProfFun;
        }
        if (std::abs(TfluxBip) < eps_heat_flux) {
          Lfac = largenum;
        }
        else {
          Lfac = -Tref_ / (kappa_ * gravity_ * TfluxBip);
        }
        double L = utau*utau*utau * Lfac;
        // limit the values of L...
        //   to be negative and finite for q>0 (unstable)
        //   to be positive and finite for q<0 (stable)
        double sgnq = (TfluxBip > 0.0)?1.0:-1.0;
        L = - sgnq * std::max(1.0e-10,std::abs(L));

        //const double theta_star = -TfluxBip / utau;
        //need TBip from temperature field - this is T(yp);
        //const double Tsurf = TBip - theta_star/kappa_*(alpha_h_*std::log(yp/z0_) + gamma_h_*yp/L);
        // evaluate rhosurf = f(Tsurf) and use rhosurf in place of rhoBip below
        double lambda = (rhoBip*kappa_*utau/(std::log(yp/z0_) - p_ABLProfFun->velocity(yp/L)));

        double tauij_nj_ui = 0.0;
        double sij_nj_ui = 0.0;
        
        double * tviscWall = stk::mesh::field_data(*tviscWall_, lfNode) ;
        double * assembledWallArea = stk::mesh::field_data(*assembledWallArea_, lfNode );
        double * assembledWallNormalDistance = stk::mesh::field_data(*assembledWallNormalDistance_, lfNode );
            
        // start the lhs assembly
        for ( int i = 0; i < nDim; ++i ) {

          double uiTan_i = 0.0;
          for ( int j = 0; j < nDim; ++j ) {
            const double ninj = p_unitNormal[i]*p_unitNormal[j];
            uiTan_i += p_uBip[i] - ninj * p_uBip[j];
            sij_nj_ui += (p_dudxBip[nDim*i+j] + p_dudxBip[nDim*j+i])*p_unitNormal[j]*p_uBip[i];            
          }
          tauij_nj_ui += lambda*(uiTan_i)*p_uBip[i];
        }

        if (sij_nj_ui < 0) {
            sij_nj_ui = 0.0;
            for (int j = 0; j < nDim; j++)
                sij_nj_ui += 0.1*p_uBip[j]*p_uBip[j]/(*assembledWallNormalDistance) ;
        }
        *tviscWall += tauij_nj_ui * aMag / (sij_nj_ui * (*assembledWallArea));        
//        *tviscWall += (sij_nj_ui * (*assembledWallArea));
        std::cerr << " Sij_nj_ui = " << sij_nj_ui <<
            "\t tauij_nj_ui = " << tauij_nj_ui <<
            "\t ttvisc = " << tauij_nj_ui * aMag / (sij_nj_ui * (*assembledWallArea)) << std::endl ;
     }
    }
  }

  std::vector<stk::mesh::FieldBase*> sum_fields(1, tviscWall_);
  stk::mesh::parallel_sum(bulk_data, sum_fields);

  for (auto part: partVec_) {
      CopyFieldAlgorithm *viscWallCopyAlg
          = new CopyFieldAlgorithm(realm_, part,
                                   tviscWall_, tvisc_,
                                   0, 1,
                                   stk::topology::NODE_RANK);
      
      viscWallCopyAlg->execute();
  }
}

} // namespace nalu
} // namespace Sierra
