/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef SinPressureAuxFunction_h
#define SinPressureAuxFunction_h

#include <AuxFunction.h>

#include <vector>

namespace sierra{
namespace nalu{

class SinPressureAuxFunction : public AuxFunction
{
public:

  SinPressureAuxFunction();

  virtual ~SinPressureAuxFunction() {}

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

class SinPressureGradientAuxFunction : public AuxFunction
{
public:

    SinPressureGradientAuxFunction(const unsigned beginPos,
                                  const unsigned endPos);

    virtual ~SinPressureGradientAuxFunction() {}

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
