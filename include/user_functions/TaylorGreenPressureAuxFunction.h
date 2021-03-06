/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef TaylorGreenPressureAuxFunction_h
#define TaylorGreenPressureAuxFunction_h

#include <AuxFunction.h>

#include <vector>

namespace sierra{
namespace nalu{

class TaylorGreenPressureAuxFunction : public AuxFunction
{
public:

  TaylorGreenPressureAuxFunction();

  virtual ~TaylorGreenPressureAuxFunction() {}
  
  virtual void do_evaluate(
    const double * coords,
    const double time,
    const unsigned spatialDimension,
    const unsigned numPoints,
    double * fieldPtr,
    const unsigned fieldSize,
    const unsigned beginPos,
    const unsigned endPos) const;
  
private:
  double uNot_;
  double pNot_;
  double rhoNot_;
  double L_;
};

} // namespace nalu
} // namespace Sierra

#endif
