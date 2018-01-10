/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeMomKEFluxEdgeOpenAlgorithm_h
#define ComputeMomKEFluxEdgeOpenAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class ComputeMomKEFluxEdgeOpenAlgorithm : public Algorithm
{
public:

  ComputeMomKEFluxEdgeOpenAlgorithm(
    Realm &realm,
    stk::mesh::Part *part);
  ~ComputeMomKEFluxEdgeOpenAlgorithm();

  void execute();

  const bool meshMotion_;

  VectorFieldType *velocity_;
  GenericFieldType *openMassFlowRate_;
};

} // namespace nalu
} // namespace Sierra

#endif
