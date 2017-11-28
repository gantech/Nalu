/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef SinVelocityAuxFunction_h
#define SinVelocityAuxFunction_h

#include <AuxFunction.h>

#include <vector>

namespace sierra{
namespace nalu{

class SinVelocityAuxFunction : public AuxFunction
{
public:

  SinVelocityAuxFunction(
    const unsigned beginPos,
    const unsigned endPos);

  virtual ~SinVelocityAuxFunction() {}

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
  const double u_infty_;
  const double alpha_;
  const double omega_;
};

} // namespace nalu
} // namespace Sierra

#endif
