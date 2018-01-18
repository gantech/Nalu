/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeMomKEFluxElemWallAlgorithm_h
#define ComputeMomKEFluxElemWallAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class ComputeMomKEFluxElemWallAlgorithm : public Algorithm
{
public:

  ComputeMomKEFluxElemWallAlgorithm(
    Realm &realm,
    stk::mesh::Part *part);
  ~ComputeMomKEFluxElemWallAlgorithm();

  void execute();

  const double includeDivU_;
  const bool meshMotion_;

  VectorFieldType *coordinates_;
  VectorFieldType *velocity_;
  ScalarFieldType *pressure_;
  ScalarFieldType *viscosity_;  
  GenericFieldType *exposedAreaVec_;

  const bool shiftMomKEFlux_;
};

} // namespace nalu
} // namespace Sierra

#endif
