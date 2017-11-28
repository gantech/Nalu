/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "user_functions/MomentumSinSrcElemKernel.h"
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
MomentumSinSrcElemKernel<AlgTraits>::MomentumSinSrcElemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  ElemDataRequests& dataPreReqs,
  bool lumped)
  : Kernel(),
    ipNodeMap_(sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_)->ipNodeMap()),
    cur_time_(0.0),
    u_infty_(10.0),
    alpha_(1.0),
    omega_(1.0),
    mu_(1.0e-5)
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
MomentumSinSrcElemKernel<AlgTraits>::~MomentumSinSrcElemKernel()
{}

template<typename AlgTraits>
void
MomentumSinSrcElemKernel<AlgTraits>::setup(const TimeIntegrator& timeIntegrator)
{
    cur_time_ = timeIntegrator.get_current_time();
}

template<typename AlgTraits>
void
MomentumSinSrcElemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType**>& lhs,
  SharedMemView<DoubleType*>& rhs,
  ScratchViews<DoubleType>& scratchViews)
{

  // Forcing nDim = 3 instead of using AlgTraits::nDim_ here to avoid compiler
  // warnings when this template is instantiated for 2-D topologies.
  DoubleType NALU_ALIGN(64) w_scvCoords[3];

  SharedMemView<DoubleType**>& v_coordinates = scratchViews.get_scratch_view_2D(*coordinates_);
  SharedMemView<DoubleType*>& v_scv_volume = scratchViews.get_me_views(CURRENT_COORDINATES).scv_volume;

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

    DoubleType x = w_scvCoords[0];
    DoubleType y = w_scvCoords[1];

    // Compute RHS contributions
    const DoubleType scV = v_scv_volume(ip);
    const int nnNdim = nearestNode * AlgTraits::nDim_;


    rhs(nnNdim + 0) += (u_infty_*((mu_ - 1.*omega_)*stk::math::cos(omega_*cur_time_ - x) - 1.*u_infty_*stk::math::sin(2*omega_*cur_time_ - 2*x) + 1.*u_infty_*stk::math::sin(2*omega_*cur_time_- 2.*x) - 1.*mu_*stk::math::sin(omega_*cur_time_ - x)) ) * scV;

    rhs(nnNdim + 1) += (mu_*u_infty_*stk::math::cos(omega_*cur_time_ - x)) * scV;

  }
}

INSTANTIATE_KERNEL(MomentumSinSrcElemKernel);

}  // nalu
}  // sierra
