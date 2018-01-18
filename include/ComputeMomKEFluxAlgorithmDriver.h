/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeMomKEFluxAlgorithmDriver_h
#define ComputeMomKEFluxAlgorithmDriver_h

#include <AlgorithmDriver.h>
#include <vector>

namespace sierra{
namespace nalu{

class Realm;
class SolutionOptions;

class ComputeMomKEFluxAlgorithmDriver : public AlgorithmDriver
{
public:

  ComputeMomKEFluxAlgorithmDriver(
    Realm &realm);

  ~ComputeMomKEFluxAlgorithmDriver();

  void compute_accumulation(std::vector<double> & , double &);
  double compute_dissipation();
  void provide_output();
  
  SolutionOptions &solnOpts_;
  bool lumpedMass_;

  void pre_work();
  void post_work();
  
};
  

} // namespace nalu
} // namespace Sierra

#endif
