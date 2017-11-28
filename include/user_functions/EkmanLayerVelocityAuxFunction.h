/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef EkmanLayerVelocityAuxFunction_h
#define EkmanLayerVelocityAuxFunction_h

#include <AuxFunction.h>

#include <vector>

namespace sierra{
namespace nalu{

class EkmanLayerVelocityAuxFunction : public AuxFunction
{
public:

  EkmanLayerVelocityAuxFunction(
    const unsigned beginPos,
    const unsigned endPos);

  virtual ~EkmanLayerVelocityAuxFunction() {}
  
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
  double G_;
  double alpha_;
  double D_;
  double pi_;

};

} // namespace nalu
} // namespace Sierra

#endif
