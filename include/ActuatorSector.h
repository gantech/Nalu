/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ActuatorSector_h
#define ActuatorSector_h

#include "Actuator.h"

namespace sierra{
namespace nalu{

class Realm;

class ActuatorSectorInfo {
public:
  ActuatorSectorInfo();
  ~ActuatorSectorInfo();

  // for each type of probe, e.g., line of site, hold some stuff
  int processorId_;
  int numPoints_;
  std::string turbineName_;
  double radius_;
  double omega_;
  double twoSigSq_;
  Coordinates tipCoordinates_;
  Coordinates tailCoordinates_;
  Coordinates coordinates_;
};

// class that holds all of the action... for each point, hold the current location and other useful info
class ActuatorSectorPointInfo {
 public:
  ActuatorSectorPointInfo(
    size_t localId, Point centroidCoords, double radius, double omega, double twoSigSq, double *velocity);
  ~ActuatorSectorPointInfo();
  size_t localId_;
  Point centroidCoords_;
  double radius_;
  double omega_;
  double twoSigSq_;
  double bestX_;
  stk::mesh::Entity bestElem_;

  // mesh motion specifics
  double velocity_[3];

  std::vector<double> isoParCoords_;
  std::vector<stk::mesh::Entity> elementVec_;
};
 
 class ActuatorSector: public Actuator
{
public:
  
  ActuatorSector(
    Realm &realm,
    const YAML::Node &node);
  ~ActuatorSector();
  
  // load all of the options
  void load(
    const YAML::Node & node);

  // setup part creation and nodal field registration (before populate_mesh())
  void setup();

  // setup part creation and nodal field registration (after populate_mesh())
  void initialize();

  // determine element bounding box in the mesh
  void populate_candidate_elements();

  // fill in the map that will hold point and ghosted elements
  void create_actuator_sector_point_info_map();

  // figure out the set of elements that belong in the custom ghosting data structure
  void determine_elems_to_ghost();

  // deal with custom ghosting
  void manage_ghosting();

  // manage rotation, now only in the y-z plane
  void set_current_coordinates(
    double *lineCentroid, double *centroidCoords, const double &omega, const double &currentTime);
  void set_current_velocity(
    double *lineCentroid, const double *centroidCoords, double *velocity, const double &omega);

  // populate vector of elements
  void complete_search();
    
  // populate nodal field and output norms (if appropriate)
  void execute();

  // support methods to gather data; scalar and vector
  void resize_std_vector( 
    const int &sizeOfField,
    std::vector<double> &theVector,   
    stk::mesh::Entity elem, 
    const stk::mesh::BulkData & bulkData);

  // general gather methods for scalar and vector (both double)
  void gather_field(
    const int &sizeOfField,
    double *fieldToFill, 
    const stk::mesh::FieldBase &stkField,
    stk::mesh::Entity const* elem_node_rels, 
    const int &nodesPerElement);

  void gather_field_for_interp(
    const int &sizeOfField,
    double *fieldToFill, 
    const stk::mesh::FieldBase &stkField,
    stk::mesh::Entity const* elem_node_rels, 
    const int &nodesPerElement);

  // element volume and scv volume populated
  double compute_volume( 
    const int &nDim,
    stk::mesh::Entity elem, 
    const stk::mesh::BulkData & bulkData);

  // interpolate field to point centroid
  void interpolate_field(
    const int &sizeOfField,
    stk::mesh::Entity elem, 
    const stk::mesh::BulkData & bulkData,
    double *isoParCoords,
    const double *fieldAtNodes,
    double *pointField);

  // drag at the point centroid
  void compute_point_drag( 
    const int &nDim, 
    const double &pointRadius,
    const double *pointVelocity,
    const double *pointGasVelocity,
    const double &pointGasViscosity,
    const double &pointGasDensity, 
    double *pointDrag,
    double &pointDragLHS);

  // centroid of the element
  void compute_elem_centroid( 
    const int &nDim,
    double *elemCentroid,
    const int &nodesPerElement);

  // radius from element centroid to point centroid
  double compute_radius( 
    const int &nDim,
    const double *elemCentroid,
    const double *pointCentroid);

  // drag fource at given radius
  void compute_elem_drag_given_radius( 
    const int &nDim,
    const double &radius, 
    const double &twoSigSq, 
    const double *pointDrag,
    double *elemDrag);

  // finally, perform the assembly
  void assemble_source_to_nodes(
    const int &nDim,
    stk::mesh::Entity elem, 
    const stk::mesh::BulkData & bulkData,
    const double &elemVolume,
    const double *drag,
    const double &dragLHS,
    stk::mesh::FieldBase &actuator_source,
    stk::mesh::FieldBase &actuator_source_lhs,
    const double &lhsFac); 

  // hold the realm
  Realm &realm_;

  // type of stk search
  stk::search::SearchMethod searchMethod_;
  
  // custom ghosting
  stk::mesh::Ghosting *actuatorSectorGhosting_;

  // how many elements to ghost?
  uint64_t needToGhostCount_;
  stk::mesh::EntityProcVec elemsToGhost_;
  
  // local id for set of points
  uint64_t localPointId_;

  // does the actuator line move?
  bool actuatorSectorMotion_;

  // everyone needs pi
  const double pi_;

  // save off product of search
  std::vector<std::pair<theKey, theKey> > searchKeyPair_;

  // bounding box data types for stk_search */
  std::vector<boundingSphere> boundingSphereVec_;
  std::vector<boundingElementBox> boundingElementBoxVec_;

  // target names for set of bounding boxes
  std::vector<std::string> searchTargetNames_;
 
  // vector of averaging information
  std::vector<ActuatorSectorInfo *> actuatorSectorInfo_; 

  // map of point info objects
  std::map<size_t, ActuatorSectorPointInfo *> actuatorSectorPointInfoMap_;

  // scratch space
  std::vector<double> ws_coordinates_;
  std::vector<double> ws_scv_volume_;
  std::vector<double> ws_velocity_;
  std::vector<double> ws_density_;
  std::vector<double> ws_viscosity_;

};


} // namespace nalu
} // namespace Sierra

#endif