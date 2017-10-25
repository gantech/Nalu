/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "Plane2D.h"
#include "FieldTypeDef.h"

#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_mesh/base/Field.hpp"
#include "stk_io/IossBridge.hpp"

namespace sierra {
namespace nalu {

Plane2D::Plane2D(
  stk::mesh::MetaData& meta,
  stk::mesh::BulkData& bulk,
  const std::string partName
) : meta_(meta),
    bulk_(bulk),
    vertices_(4, std::vector<double>(3, 0.0)),
    blockName_(partName)
{}

void Plane2D::setup()
{
  stk::mesh::Part* part = meta_.get_part(blockName_);

  if (part != nullptr) {
    throw std::runtime_error(
      "Plane2D::setup(): Cannot overwrite existing part in mesh database: " + blockName_);
  }
  else {
    stk::mesh::Part& newPart = meta_.declare_part(blockName_, stk::topology::NODE_RANK);
    stk::io::put_io_part_attribute(newPart);
  }
}

void Plane2D::initialize()
{
  const size_t iproc = bulk_.parallel_rank();
  const size_t nproc = bulk_.parallel_size();
  stk::mesh::Part& part = *meta_.get_part(blockName_);

  size_t gNumPoints = nx_ * ny_;
  unsigned numPoints = 0;
  unsigned offset = 0;

  if ((gNumPoints < nproc) && (iproc < gNumPoints)) {
    numPoints = 1;
  } else {
    numPoints = gNumPoints / nproc;
    offset = iproc * numPoints;
    unsigned rem = gNumPoints % nproc;

    if ((rem > 0) && (iproc < rem)) numPoints++;
    offset += (iproc < rem)? iproc : rem;
  }

  std::vector<stk::mesh::EntityId> newIDs(numPoints);
  std::vector<stk::mesh::Entity> nodeVec(numPoints);

  bulk_.modification_begin();
  if (numPoints > 0) {
    bulk_.generate_new_ids(stk::topology::NODE_RANK, numPoints, newIDs);

    for(unsigned i=0; i<numPoints; i++) {
      stk::mesh::Entity node = bulk_.declare_entity(
        stk::topology::NODE_RANK, newIDs[i], part);
      nodeVec[i] = node;
    }
  }
  bulk_.modification_end();

  VectorFieldType* coords = meta_.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, "coordinates");

  for (unsigned k=0; k < numPoints; k++) {
    int j = (offset+k) / nx_;
    int i = (offset+k) % nx_;
    double* pt = stk::mesh::field_data(*coords, nodeVec[k]);

    const double rx = i * dx_;
    const double ry = j * dy_;

    pt[0] = ((1.0 - rx) * (1.0 - ry) * vertices_[0][0] +
             rx * (1.0 - ry) * vertices_[1][0] +
             rx * ry * vertices_[2][0] +
             (1.0 - rx) * ry * vertices_[3][0]);
    pt[1] = ((1.0 - rx) * (1.0 - ry) * vertices_[0][1] +
             rx * (1.0 - ry) * vertices_[1][1] +
             rx * ry * vertices_[2][1] +
             (1.0 - rx) * ry * vertices_[3][1]);
    pt[2] = ((1.0 - rx) * (1.0 - ry) * vertices_[0][2] +
             rx * (1.0 - ry) * vertices_[1][2] +
             rx * ry * vertices_[2][2] +
             (1.0 - rx) * ry * vertices_[3][2]);
  }
}

}  // nalu
}  // sierra
