/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeMomKEFluxElemOpenAlgorithm_h
#define ComputeMomKEFluxElemOpenAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class ComputeMomKEFluxElemOpenAlgorithm : public Algorithm
{
public:

  ComputeMomKEFluxElemOpenAlgorithm(
    Realm &realm,
    stk::mesh::Part *part);
  ~ComputeMomKEFluxElemOpenAlgorithm();

  void execute();

  const bool meshMotion_;

  VectorFieldType *velocity_;
  GenericFieldType *exposedAreaVec_;
  GenericFieldType *openMassFlowRate_;

  const bool shiftMomKEFlux_;
};

} // namespace nalu
} // namespace Sierra

#endif
