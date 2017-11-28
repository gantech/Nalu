/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/SinPressureAuxFunction.h>
#include <algorithm>

// basic c++
#include <cmath>
#include <vector>
#include <stdexcept>

namespace sierra{
namespace nalu{

SinPressureAuxFunction::SinPressureAuxFunction() :
  AuxFunction(0,1),
  u_infty_(10.0),
  alpha_(1.0),
  omega_(1.0)
{
  // does nothing
}

void
SinPressureAuxFunction::do_evaluate(
  const double *coords,
  const double time,
  const unsigned spatialDimension,
  const unsigned numPoints,
  double * fieldPtr,
  const unsigned fieldSize,
  const unsigned /*beginPos*/,
  const unsigned /*endPos*/) const
{

  for(unsigned p=0; p < numPoints; ++p) {

    const double x = coords[0];
    const double y = coords[1];

    fieldPtr[0] = u_infty_*u_infty_*std::cos(omega_*time - x)*std::cos(omega_*time - x)* + 2.*u_infty_*u_infty_*std::cos(x)*std::sin(omega_*time) - 2.0*u_infty_*u_infty_*std::cos(omega_*time)*std::sin(x);

    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

SinPressureGradientAuxFunction::SinPressureGradientAuxFunction(
   const unsigned beginPos,
   const unsigned endPos) :
   AuxFunction(beginPos, endPos),
    u_infty_(10.0),
    alpha_(1.0),
    omega_(1.0)
{
    // does nothing
}
    
void
SinPressureGradientAuxFunction::do_evaluate(
  const double *coords,
  const double time,
  const unsigned spatialDimension,
  const unsigned numPoints,
  double * fieldPtr,
  const unsigned fieldSize,
  const unsigned /*beginPos*/,
  const unsigned /*endPos*/) const
{
    
    for(unsigned p=0; p < numPoints; ++p) {
        
        const double x = coords[0];
        const double y = coords[1];

        fieldPtr[0] = -2.0* u_infty_*u_infty_* std::cos(omega_* time - x) * (1 - std::sin(omega_* time - x));
        fieldPtr[1] = 0.0;
        
        fieldPtr += fieldSize;
        coords += spatialDimension;
    }
}
    
} // namespace nalu
} // namespace Sierra
