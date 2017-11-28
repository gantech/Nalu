/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/EkmanLayerVelocityAuxFunction.h>
#include <algorithm>

// basic c++
#include <cmath>
#include <vector>
#include <stdexcept>

namespace sierra{
namespace nalu{

EkmanLayerVelocityAuxFunction::EkmanLayerVelocityAuxFunction(
  const unsigned beginPos,
  const unsigned endPos) :
  AuxFunction(beginPos, endPos),
  G_(1.0), //15.0
  alpha_(0.0),
  D_(1.0), //261.85340212559663
  pi_(std::acos(-1.0))
{
  // nothing to do
}

void
EkmanLayerVelocityAuxFunction::do_evaluate(
  const double *coords,
  const double t,
  const unsigned /*spatialDimension*/,
  const unsigned numPoints,
  double * fieldPtr,
  const unsigned fieldSize,
  const unsigned /*beginPos*/,
  const unsigned /*endPos*/) const
{
  for(unsigned p=0; p < numPoints; ++p) {

    double x = coords[0];
    double y = coords[1];
    double z = coords[2];

    fieldPtr[0] = G_ * ( cos(alpha_) - exp(-z/D_)*cos(z/D_ - alpha_) );
    fieldPtr[1] = G_ * ( sin(alpha_) + exp(-z/D_)*sin(z/D_ - alpha_) );
    fieldPtr[2] = 0.0;

    fieldPtr += fieldSize;
    coords += fieldSize;
  }
}

} // namespace nalu
} // namespace Sierra
