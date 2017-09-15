/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef PLANE2D_H
#define PLANE2D_H

#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/BulkData.hpp"

#include "yaml-cpp/yaml.h"

#include <string>

namespace sierra {
namespace nalu {

/** Planar 2D nodeset generation utility
 *
 *  Given a part name, lateral and longitudinal dimensions, and a set of
 *  vertices indicating the corners of the 2-D quadrilateral plane, this class
 *  generates a spatial set of grid points that can be used to sample data off
 *  of the computational mesh. This class is meant to be used from within other
 *  Nalu classes such as sierra::nalu::DataProbePostProcessing
 *
 *  \sa ABLForcingAlgorithm
 */
class Plane2D
{
public:
  /**
   *  @param meta STK MetaData instance
   *  @param bulk STK BulkData instance
   *  @param partName Name of the nodeset to be generated
   *
   */
  Plane2D(
    stk::mesh::MetaData& meta,
    stk::mesh::BulkData& bulk,
    const std::string partName);

  ~Plane2D() {}

  //! Register parts to STK meta data object
  void setup();

  //! Create nodal entities and populate the coordinates
  void initialize();

  /** Get/set the internal vertices describing the 2-D plane
   *
   *  This is an array of shape (4, 3) that describes the extents of the 2-D
   *  plane in 3-D space. Users are expected to populate the vertices
   *  appropriately before calling the Plane2D::initialize() method.
   */
  inline std::vector<std::vector<double>>& vertices()
  { return vertices_; }

  /** Set the dimensions used to discretize the plane into a grid of nodal points
   *
   *  @param mx[in] Dimensions in the `x` direction
   *  @param my[in] Dimension in the `y` direction
   */
  inline void set_dimensions(size_t mx, size_t my)
  {
    mx_ = mx;
    my_ = my;
    nx_ = mx_ + 1;
    ny_ = my_ + 1;

    dx_ = 1.0 / static_cast<double>(mx_);
    dy_ = 1.0 / static_cast<double>(my_);
  }

private:
  Plane2D() = delete;
  Plane2D(const Plane2D&) = delete;

  //! STK Metadata object
  stk::mesh::MetaData& meta_;

  //! STK Bulkdata object
  stk::mesh::BulkData& bulk_;

  //! Corners of the computational domain mesh
  std::vector<std::vector<double>> vertices_;

  //! Name of the part or nodeset
  std::string blockName_{"block_1"};

  //! Spatial resolutions
  double dx_, dy_;

  //! Mesh sizes
  size_t nx_, ny_;
  size_t mx_, my_;
};

}  // nalu
}  // sierra

#endif /* PLANE2D_H */
