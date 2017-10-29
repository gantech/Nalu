/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifdef NALU_USES_SMD

/** @file ActuatorLineSMD.h
 *  @brief A class to couple Nalu with SMD for actuator line simulation of a spring-mass-damper system
 *
 */

#ifndef ActuatorLineSMD_h
#define ActuatorLineSMD_h

#include <stk_util/parallel/ParallelVectorConcat.hpp>
#include "Actuator.h"

// SMD C++ API
#include "smd.H"

namespace sierra{
namespace nalu{

class Realm;

/** Class that holds all of the information relevant to each spring-mass-damper system
 *
 *
 */
//
class ActuatorLineSMDInfo {
public:
  ActuatorLineSMDInfo();
  ~ActuatorLineSMDInfo();

  int processorId_; ///< The processor on which SMD is run for this spring-mass-damper system
  int numPoints_; ///< The total number of actuator points for this spring-mass-damper system
  std::string smdName_; ///< The spring-mass-damper system name
  Coordinates epsilon_; ///< The Gaussian spreading width in (chordwise, spanwise, thickness) directions
};

/** Class that holds all of the search action for each actuator point
 *
 *
 */
//
class ActuatorLineSMDPointInfo {
 public:
  ActuatorLineSMDPointInfo(
			    size_t globSMDId, Point centroidCoords, double searchRadius, Coordinates epsilon);
  ~ActuatorLineSMDPointInfo();
  size_t globSMDId_; ///< Global spring-mass-damper system number.
  Point centroidCoords_; ///< The coordinates of the actuator point.
  double searchRadius_; ///< Elements within this search radius will be affected by this actuator point.
  Coordinates epsilon_; ///< The Gaussian spreading width in (chordwise, spanwise, thickness) directions for this actuator point.
  double bestX_; ///< A number returned by stk::isInElement that determines whether an actuator point is inside (< 1) or outside an element (> 1). However, we choose the bestElem_ for this actuator point to be the one with the lowest bestX_.
  stk::mesh::Entity bestElem_; ///< The element within which the actuator point lies.

  std::vector<double> isoParCoords_; ///< The isoparametric coordinates of the bestElem_.
  std::vector<stk::mesh::Entity> elementVec_; ///< A list of elements that lie within the searchRadius_ around the actuator point.
};

/** The ActuatorLineSMD class couples Nalu with the third party library SMD for actuator line simulations of a spring-mass-damper system
 *
 * SMD is a library that models a point spring-mass-damper system. The effect of the turbine
 * on the flow field is modeled using the actuator line approach. The force exerted by the wind
 * turbine on the flow field is lumpled into a set of body forces at a discrete set of actuator
 * points. This class spreads the the body force at each actuator point using a Gaussian function.

 * 1) During the load phase - the turbine data from the yaml file is read and stored in an
 *    object of the smd::smdInputs class

 * 2) During the initialize phase - The processor containing the hub of each turbine is found
 *    through a search and assigned to be the one controlling SMD for that turbine. All
 *    processors controlling > 0 turbines initialize FAST, populate the map of ActuatorLinePointInfo
 *    and initialize element searches for all the actuator points associated with the turbines.
 *
 * 3) Elements are ghosted to the owning point rank. We tried the opposite approach of
 *    ghosting the actuator points to the processor owning the elements. The second approach
 *    was found to peform poorly compared to the first method.
 *
 * 4) A time lagged simple FSI model is used to interface Nalu with the turbine model:
 *    + The velocity at time step at time step 'n' is sampled at the actuator points and sent
 *       to SMD
 *    + SMD advances the spring-mass-damper system upto the next Nalu time step 'n+1'
 *    + The body forces at the actuator points are converted to the source terms of the momentum
 *      equation to advance Nalu to the next time step 'n+1'.
 *
 * 5) During the execute phase called every time step, we sample the velocity at the spring-mass-damper
 *    point and pass it to SMD. The spring-mass-damper system is advanced upto Nalu's
 *    next time step to get the body forces at the actuator point. We then iterate over the
 *    ActuatorLinePointInfoMap (now only single point) to assemble source terms.
 *
 *    actuator:
 *     type: ActLineSMD
 *     search_method: boost_rtree
 *     search_target_part: Unspecified-2-HEX
 *     dry_run:  False
 *     debug:    False
 *     t_max:    5.0
 *     n_every_checkpoint: 100
 *     epsilon: [ 5.0, 5.0, 5.0 ]
 */

class ActuatorLineSMD: public Actuator {
 public:

  ActuatorLineSMD(
    Realm &realm,
    const YAML::Node &node);
  ~ActuatorLineSMD();

  // load all of the options
  void load(
    const YAML::Node & node);

  // load the options for each turbine
  void readSMDData(int iTurb, smd::smdInputs & fi, YAML::Node turbNode);

  // setup part creation and nodal field registration (before populate_mesh())
  void setup();

  // allocate smd to processor containing base location
  void allocateSMDToProc() ;

  // Allocate SMD to to processor0, initialize SMD and get location of actuator points
  void initialize();

  // setup part creation and nodal field registration (after populate_mesh())
  void update();

  // determine processor bounding box in the mesh
  void populate_candidate_procs();

  // determine element bounding box in the mesh
  void populate_candidate_elements();

  // fill in the map that will hold point and ghosted elements
  void create_actuator_line_point_info_map();

  // figure out the set of elements that belong in the custom ghosting data structure
  void determine_elems_to_ghost();

  // deal with custom ghosting
  void manage_ghosting();

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

  void gather_field_inverted_order(
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
    const std::vector<double> & elemForce,
    const double &gLocal,
    stk::mesh::FieldBase & elemCoords,
    stk::mesh::FieldBase &actuator_source,
    stk::mesh::FieldBase &g,
    stk::mesh::FieldBase &dualNodalVolume,
    std::vector<double>& forceSum);

  Realm &realm_; ///< hold the realm

  stk::search::SearchMethod searchMethod_; ///< type of stk search

  stk::mesh::Ghosting *actuatorLineGhosting_;  ///< custom ghosting
  uint64_t needToGhostCount_;  ///< how many elements to ghost?
  stk::mesh::EntityProcVec elemsToGhost_; ///< elements to ghost

  int tStepRatio_;  ///< Ratio of Nalu time step to SMD time step (dtNalu/dtSMD) - Should be an integral number

  std::vector<std::pair<theKey, theKey> > searchKeyPair_;  ///< save off product of search

  // bounding box data types for stk_search
  std::vector<boundingSphere> boundingSphereVec_; ///< bounding box around each actuator point
  std::vector<boundingElementBox> boundingElementBoxVec_; ///< bounding box around elements
  std::vector<boundingSphere> boundingHubSphereVec_; ///< bounding box around the hub point of each turbine
  std::vector<boundingElementBox> boundingProcBoxVec_; ///< bounding box around all the nodes residing locally on each processor

  std::vector<std::string> searchTargetNames_;  ///< target names for set of bounding boxes

  std::vector<ActuatorLineSMDInfo *> actuatorLineInfo_;   ///< vector of objects containing information for each turbine

  std::map<size_t, ActuatorLineSMDPointInfo *> actuatorLinePointInfoMap_;  ///< map of point info objects

  // scratch space
  std::vector<double> ws_coordinates_;
  std::vector<double> ws_scv_volume_;
  std::vector<double> ws_velocity_;
  std::vector<double> ws_density_;
  std::vector<double> ws_viscosity_;

  smd::smdInputs i_smd; ///< Object to hold input information for SMD
  smd::smd p_smd; ///< SMD API handle

  std::vector<std::vector<double>> thrust;
  std::vector<std::vector<double>> torque;

};


} // namespace nalu
} // namespace Sierra

#endif

#endif
