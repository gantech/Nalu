/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef CONTINUITYSMDSRCELEMKERNEL_H
#define CONTINUITYSMDSRCELEMKERNEL_H

#include "Kernel.h"
#include "FieldTypeDef.h"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <Kokkos_Core.hpp>

namespace sierra {
namespace nalu {

class SolutionOptions;
class MasterElement;
class ElemDataRequests;

/** CMM source term for MMS test for Spring-Mass-Damper system using actuator lines for continuty equation (pressure DOF)
 */
template<typename AlgTraits>
class ContinuitySMDSrcElemKernel: public Kernel
{
public:
  ContinuitySMDSrcElemKernel(
    const stk::mesh::BulkData&,
    SolutionOptions&,
    ElemDataRequests&);

  virtual ~ContinuitySMDSrcElemKernel() {}

  /** Perform pre-timestep work for the computational kernel
   */
  virtual void setup(const TimeIntegrator&);

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<DoubleType**>&,
    SharedMemView<DoubleType*>&,
    ScratchViews<DoubleType>&);

private:
  ContinuitySMDSrcElemKernel() = delete;

  VectorFieldType *coordinates_;

  const int *ipNodeMap_;

  DoubleType cur_time_;
  const DoubleType u_infty_;
  const DoubleType A_;
  const DoubleType sigma_;
  const DoubleType alpha_;
  const DoubleType omega_;

  // fixed scratch space
  Kokkos::View<DoubleType[AlgTraits::numScvIp_][AlgTraits::nodesPerElement_]> v_shape_function_{"v_shape_function"};
};

}  // nalu
}  // sierra

#endif /* CONTINUITYSMDSRCELEMKERNEL_H */
