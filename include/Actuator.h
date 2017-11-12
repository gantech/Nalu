/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef Actuator_h
#define Actuator_h

#include <NaluParsing.h>
#include<FieldTypeDef.h>

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Ghosting.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_search/SearchMethod.hpp>

// stk forwards
/*namespace stk {
  namespace mesh {
    struct Entity;
  }
  }*/

// basic c++
#include <string>
#include <vector>
#include <utility>

namespace sierra{
namespace nalu{

// common type defs
typedef stk::search::IdentProc<uint64_t,int> theKey;
typedef stk::search::Point<double> Point;
typedef stk::search::Sphere<double> Sphere;
typedef stk::search::Box<double> Box;
typedef std::pair<Sphere,theKey> boundingSphere;
typedef std::pair<Box,theKey> boundingElementBox;

class Realm;
 
class Actuator
{
public:
  
  Actuator(
    Realm &realm,
    const YAML::Node &node) {} 
  virtual ~Actuator() {}
  
  // load all of the options
  virtual void load(
    const YAML::Node & node) = 0;

  // setup part creation and nodal field registration (before populate_mesh())
  virtual void setup() = 0;

  // setup part creation and nodal field registration (after populate_mesh())
  virtual void initialize() = 0;

  // predict the state of the structural model at the next time step
  virtual void predict_struct_time_step() = 0;

  // firmly advance the state of the structural model to the next time step
  virtual void advance_struct_time_step() = 0;
  
  // sample velocity at the actuator points and send to the structural model
  virtual void sample_vel() = 0;

  // populate nodal field and output norms (if appropriate)
  virtual void execute() = 0;

};

} // namespace nalu
} // namespace Sierra

#endif
