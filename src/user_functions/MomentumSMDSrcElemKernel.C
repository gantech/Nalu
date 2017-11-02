/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "user_functions/MomentumSMDSrcElemKernel.h"
#include "AlgTraits.h"
#include "master_element/MasterElement.h"
#include "SolutionOptions.h"
#include "TimeIntegrator.h"

// template and scratch space
#include "BuildTemplates.h"
#include "ScratchViews.h"

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

namespace sierra {
namespace nalu {

template<typename AlgTraits>
MomentumSMDSrcElemKernel<AlgTraits>::MomentumSMDSrcElemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  ElemDataRequests& dataPreReqs,
  bool lumped)
  : Kernel(),
    ipNodeMap_(sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_)->ipNodeMap()),
    cur_time_(0.0),
    u_infty_(2.0),
    A_(1.0),
    sigma_(0.5),
    alpha_(1.0),
    omega_(1.0),
    M_(10.0),
    C_(1.0),
    K_(1.0),
    mu_(1e-5)
{
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  coordinates_ = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, solnOpts.get_coordinates_name());

  MasterElement* meSCV = sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_);

  if ( lumped ) {
    get_scv_shape_fn_data<AlgTraits>([&](double* ptr){meSCV->shifted_shape_fcn(ptr);}, v_shape_function_);
  }
  else {
    get_scv_shape_fn_data<AlgTraits>([&](double* ptr){meSCV->shape_fcn(ptr);}, v_shape_function_);
  }

  // add master elements
  dataPreReqs.add_cvfem_volume_me(meSCV);

  // fields and data
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);
  dataPreReqs.add_master_element_call(SCV_VOLUME, CURRENT_COORDINATES);
}

template<typename AlgTraits>
MomentumSMDSrcElemKernel<AlgTraits>::~MomentumSMDSrcElemKernel()
{}

template<typename AlgTraits>
void
MomentumSMDSrcElemKernel<AlgTraits>::setup(const TimeIntegrator& timeIntegrator)
{
    cur_time_ = timeIntegrator.get_current_time();
}

template<typename AlgTraits>
void
MomentumSMDSrcElemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType**>& lhs,
  SharedMemView<DoubleType*>& rhs,
  ScratchViews<DoubleType>& scratchViews)
{

  // Forcing nDim = 3 instead of using AlgTraits::nDim_ here to avoid compiler
  // warnings when this template is instantiated for 2-D topologies.
  DoubleType w_scvCoords[3] KOKKOS_ALIGN(64);

  SharedMemView<DoubleType**>& v_coordinates = scratchViews.get_scratch_view_2D(*coordinates_);
  SharedMemView<DoubleType*>& v_scv_volume = scratchViews.get_me_views(CURRENT_COORDINATES).scv_volume;

  DoubleType sigma2 = sigma_ * sigma_;
  DoubleType oneOverSigma2 = 1.0/sigma2 ;
  DoubleType oneOverSigma4 = oneOverSigma2 * oneOverSigma2 ;

  for (int ip=0; ip < AlgTraits::numScvIp_; ++ip) {

    // nearest node to ip
    const int nearestNode = ipNodeMap_[ip];

    // zero out
    for ( int j =0; j < AlgTraits::nDim_; ++j )
        w_scvCoords[j] = 0.0;

    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
        const DoubleType r = v_shape_function_(ip,ic);
        for ( int j = 0; j < AlgTraits::nDim_; ++j )
            w_scvCoords[j] += r*v_coordinates(ic,j);
    }

    DoubleType xCoord = w_scvCoords[0];
    DoubleType yCoord = w_scvCoords[1];

    DoubleType yrels = yCoord - alpha_ * stk::math::sin(omega_ * cur_time_);
    DoubleType expfn = stk::math::exp(-(xCoord*xCoord + yrels*yrels)*oneOverSigma2);

    // Compute RHS contributions
    const DoubleType scV = v_scv_volume(ip);
    const int nnNdim = nearestNode * AlgTraits::nDim_;
    rhs(nnNdim + 0) += (A_*expfn*oneOverSigma2*(-2*A_*xCoord + 3.*(-(A_*expfn) + u_infty_)*xCoord - 2.*alpha_*omega_*stk::math::cos(omega_*cur_time_)*yrels - (2.*A_*expfn*xCoord*yrels*yrels)*oneOverSigma2 + (2.*(-(A_*expfn) + u_infty_)*xCoord*yrels*yrels)*oneOverSigma2 - mu_*(4.0 - (4.0*xCoord*xCoord)*oneOverSigma2 - (4*yrels*yrels)*oneOverSigma2 - (2*(xCoord + (2*xCoord*yrels*yrels)*oneOverSigma2))/3.)))*scV ;

    rhs(nnNdim + 1) += (A_*expfn*oneOverSigma2*(1.*alpha_*omega_*xCoord*stk::math::cos(omega_*cur_time_) - 2*A_*yrels - 1.*(-(A_*expfn) + u_infty_)*yrels + (2.*(-(A_*expfn) + u_infty_)*xCoord*xCoord*yrels)*oneOverSigma2 - (2.*alpha_*omega_*xCoord*stk::math::cos(omega_*cur_time_)*yrels*yrels)*oneOverSigma2 - (4.*A_*expfn*xCoord*xCoord*yrels*yrels*yrels)*oneOverSigma4 - (alpha_*sigma2*(C_*omega_*stk::math::cos(omega_*cur_time_) + (K_ - M_*omega_*omega_)*stk::math::sin(omega_*cur_time_)))/A_ - mu_*((12.0*xCoord*yrels)*oneOverSigma2 - (4.0*xCoord*xCoord*xCoord*yrels)*oneOverSigma4 - (4.0*xCoord*yrels*yrels*yrels)*oneOverSigma4 - (2.0*(xCoord + (2.0*xCoord*yrels*yrels)*oneOverSigma2))/3.)))*scV;

  }
}

INSTANTIATE_KERNEL(MomentumSMDSrcElemKernel);

}  // nalu
}  // sierra
