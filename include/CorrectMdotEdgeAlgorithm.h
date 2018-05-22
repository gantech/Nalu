/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef CorrectMdotEdgeAlgorithm_h
#define CorrectMdotEdgeAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class CorrectMdotEdgeAlgorithm : public Algorithm
{
public:

  CorrectMdotEdgeAlgorithm(
    Realm &realm,
    stk::mesh::Part *part);
  ~CorrectMdotEdgeAlgorithm();

  void execute();
  
  const bool meshMotion_;
  VectorFieldType *velocityRTM_;
  VectorFieldType *Gpdx_;
  VectorFieldType *coordinates_;
  ScalarFieldType *pressure_;
  ScalarFieldType *density_;
  VectorFieldType *edgeAreaVec_;
  ScalarFieldType *massFlowRate_;

};

} // namespace nalu
} // namespace Sierra

#endif
