/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeMomKEFluxInflowAlgorithm_h
#define ComputeMomKEFluxInflowAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class ComputeMomKEFluxInflowAlgorithm : public Algorithm
{
public:

  ComputeMomKEFluxInflowAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    bool useShifted);
  ~ComputeMomKEFluxInflowAlgorithm();

  void execute();
  
  const double includeDivU_;
  const bool useShifted_;

  VectorFieldType *coordinates_;
  VectorFieldType *velocity_;
  ScalarFieldType *viscosity_;  
  ScalarFieldType *density_;
  ScalarFieldType *pressure_;
  GenericFieldType *exposedAreaVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
