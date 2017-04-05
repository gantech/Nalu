/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ActuatorLineFAST_h
#define ActuatorLineFAST_h

#include <stk_util/parallel/ParallelVectorConcat.hpp>
#include "Actuator.h"

// FAST c interface
#include "FAST_cInterface.h"



namespace sierra{
namespace nalu{

class Realm;

struct minDist {
  double value;
  int id;
};

// class that holds all of the search action... for each actuator point
class ActuatorLineFASTPointInfo {
 public:
  ActuatorLineFASTPointInfo(double searchRadius);
  ~ActuatorLineFASTPointInfo();
  double searchRadius_;
  stk::mesh::Entity bestElem_;
  std::vector<double> isoParCoords_;
  std::vector<stk::mesh::Entity> elementVec_;
};

class ActuatorLineFASTInfo {
public:
  ActuatorLineFASTInfo();
  ~ActuatorLineFASTInfo();

  // for each type of probe, e.g., line of site, hold some stuff
  int procId_;  // Id of the processor owning FAST in global comm
  size_t globTurbId_; // Global turbine number
  std::string turbineName_;

  int procIdGroup_; // Id of the processor owning FAST in turbine group comm
  MPI_Group turbineGroup_;
  MPI_Comm turbineComm_;

  int numPoints_;
  Coordinates epsilon_;

  double ** coords_; // Coordinates of the actuator points
  double ** velocity_; // Sampled velocity at the actuator points
  double ** force_; // Body force at the actuator points
  int * nType_; // Actuator node type of each actuator point

  minDist * bestX_; // Distance of closest element to each actuator point
  
  std::map<size_t, ActuatorLineFASTPointInfo *> actuatorLinePointInfoMap_; // Map of point info objects
  // bounding box data types for stk_search */
  std::vector<boundingSphere> boundingSphereVec_;

};

class ActuatorLineFAST: public Actuator
{
public:
  
  ActuatorLineFAST(
    Realm &realm,
    const YAML::Node &node);
  ~ActuatorLineFAST();
  
  // load all of the options
  void load(
    const YAML::Node & node);

  // setup part creation and nodal field registration (before populate_mesh())
  void setup();

  // allocate turbines to processors containing hub location
  void allocateTurbinesToProcs() ;
  
  // Allocate turbines to processors, initialize FAST and get location of actuator points
  void initialize();

  // setup part creation and nodal field registration (after populate_mesh())
  void update();

  // determine processor bounding box in the mesh
  void populate_candidate_procs();

  // determine element bounding box in the mesh
  void populate_candidate_elements();

  // fill in the map that will hold point and ghosted elements
  void create_actuator_line_point_info_map();

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

  // centroid of the element
  void compute_elem_centroid( 
    const int &nDim,
    double *elemCentroid,
    const int &nodesPerElement);

  // distance from element centroid to point centroid
  double compute_distance( 
    const int &nDim,
    const double *elemCentroid,
    const double *pointCentroid);

  // compute the body force at an element given a
  // projection weighting.
  void compute_elem_force_given_weight(
    const int &nDim,
    const double &g,
    const double *pointForce,
    double *elemForce);

  // isotropic Gaussian projection function.
  double isotropic_Gaussian_projection(
    const int &nDim,
    const double &dis,
    const Coordinates &epsilon);

  // finally, perform the assembly
  void assemble_source_to_nodes(
    const int &nDim,
    stk::mesh::Entity elem,
    const stk::mesh::BulkData & bulkData,
    const double &elemVolume,
    const double *drag,
    const double &dragLHS,
    const double &gLocal,
    stk::mesh::FieldBase &actuator_source,
    stk::mesh::FieldBase &actuator_source_lhs,
    stk::mesh::FieldBase &g,
    const double &lhsFac);

  // hold the realm
  Realm &realm_;

  // type of stk search
  stk::search::SearchMethod searchMethod_;
  
  // does the actuator line move?
  bool actuatorLineMotion_;

  // save off product of search
  std::vector<std::pair<theKey, theKey> > searchKeyPair_;

  // bounding box data types for stk_search */
  std::vector<boundingElementBox> boundingElementBoxVec_;
  std::vector<boundingSphere> boundingHubSphereVec_;
  std::vector<boundingSphere> boundingHubBigSphereVec_;
  std::vector<boundingElementBox> boundingProcBoxVec_;

  // target names for set of bounding boxes
  std::vector<std::string> searchTargetNames_;
 
  // vector of ActuatorLineFASTInfo objects
  std::vector<ActuatorLineFASTInfo *> actuatorLineInfo_; 

  // scratch space
  std::vector<double> ws_coordinates_;
  std::vector<double> ws_scv_volume_;
  std::vector<double> ws_velocity_;
  std::vector<double> ws_density_;
  std::vector<double> ws_viscosity_;

  // FAST cInterface handle
  FAST_cInterface FAST;
  
  MPI_Group worldMPIGroup_;

};


} // namespace nalu
} // namespace Sierra

#endif
