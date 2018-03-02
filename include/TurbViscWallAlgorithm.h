/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef TurbViscWallAlgorithm_h
#define TurbViscWallAlgorithm_h

#include<Algorithm.h>

#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class TurbViscWallAlgorithm : public Algorithm
{
public:
  
  TurbViscWallAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    const bool &useShifted,
    const double &gravity,
    const double &z0,
    const double &Tref);
  virtual ~TurbViscWallAlgorithm() {}
  virtual void execute();

  const bool useShifted_;
  const double z0_;
  const double Tref_;
  const double gravity_;
  const double alpha_h_;
  const double beta_m_;
  const double beta_h_;
  const double gamma_m_;
  const double gamma_h_;
  const double kappa_;
  
  ScalarFieldType *bcHeatFlux_;
  ScalarFieldType *specificHeat_;
  ScalarFieldType *density_;
  VectorFieldType *velocity_;  
  ScalarFieldType *tvisc_;
  ScalarFieldType *tviscWall_;
  ScalarFieldType *assembledWallArea_;
  GenericFieldType *wallFrictionVelocityBip_;
  GenericFieldType *wallNormalDistanceBip_;
  GenericFieldType *dudx_;
  GenericFieldType *exposedAreaVec_;  
  ScalarFieldType *assembledWallNormalDistance_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
