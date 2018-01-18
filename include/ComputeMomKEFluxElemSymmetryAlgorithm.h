/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeMomKEFluxElemSymmetryAlgorithm_h
#define ComputeMomKEFluxElemSymmetryAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class ComputeMomKEFluxElemSymmetryAlgorithm : public Algorithm
{
public:

  ComputeMomKEFluxElemSymmetryAlgorithm(
    Realm &realm,
    stk::mesh::Part *part);
  ~ComputeMomKEFluxElemSymmetryAlgorithm();

  void execute();

  const double includeDivU_;

  VectorFieldType *coordinates_;
  VectorFieldType *velocity_;
  ScalarFieldType *viscosity_;
  ScalarFieldType *pressure_;    
  GenericFieldType *exposedAreaVec_;

};

} // namespace nalu
} // namespace Sierra

#endif
