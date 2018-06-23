/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef CorrectMdotEdgeOpenAlgorithm_h
#define CorrectMdotEdgeOpenAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class CorrectMdotEdgeOpenAlgorithm : public Algorithm
{
public:

  CorrectMdotEdgeOpenAlgorithm(
    Realm &realm,
    stk::mesh::Part *part);
  ~CorrectMdotEdgeOpenAlgorithm();

  void execute();

  const bool meshMotion_;

  VectorFieldType *velocityRTM_;
  VectorFieldType *uDiagInv_;
  VectorFieldType *Gpdx_;
  VectorFieldType *coordinates_;
  ScalarFieldType *pressure_;
  ScalarFieldType *density_;
  GenericFieldType *exposedAreaVec_;
  GenericFieldType *openMassFlowRate_;
  ScalarFieldType *pressureBc_;
};

} // namespace nalu
} // namespace Sierra

#endif
