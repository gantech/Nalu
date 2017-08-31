
#ifndef ABLPOSTPROCESSINGALGORITHM_H
#define ABLPOSTPROCESSINGALGORITHM_H

#include "NaluParsing.h"
#include "FieldTypeDef.h"

#include "stk_mesh/base/Selector.hpp"

#include <string>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <unordered_set>

namespace sierra {
namespace nalu {

class Realm;
class Transfer;
class Transfers;
class SpatialAveragingAlgorithm;
/**
 * \brief ABL Forcing Source terms for Momentum and Temperature equations
 *
 * This class parses the user inputs and provides an planar statistics based 
 * postprocessing implementation within Nalu.
 * The ABL postprocessing capability is turned on by the presence of a sub-section
 * titled `abl_postprocessing` within the Realm section of the Nalu input file.
 *
 * ```
 *   abl_postprocessing:
 *     search_method: stk_kdtree
 *     search_tolerance: 0.0001
 *     search_expansion_factor: 1.5
 *     from_target_part: [Unspecified-2-HEX]
 *     target_part_format: "zplane_%.1f"
 *     heights: [80.0]
 * ```
 *
 */
class ABLPostProcessingAlgorithm
{
public:
  template <typename T>
  using Array2D = std::vector<std::vector<T>>;

  ABLPostProcessingAlgorithm(Realm&, const YAML::Node&);

  ABLPostProcessingAlgorithm(Realm&, const YAML::Node&, SpatialAveragingAlgorithm& spatialAvg);
  
  ~ABLPostProcessingAlgorithm();

  //! Parse input file for user options and initialize
  void load(const YAML::Node&);

  //! Setup ABL postprocessing (steps before mesh creation)
  void setup();

  //! Initialize ABL postprocessing (steps after mesh creation)
  void initialize();

  //! Execute field transfers, compute planar averaging, and determine source
  //! terms at desired levels.
  void execute();

  //! Evaluate the ABL postprocessing source contribution at a node
  void eval_vel_mean(
    const double,        //!< Height of the node from terrain
    std::vector<double>& //!< Source vector to be populated
    );

  //! Evaluate the ABL postprocessing source contribution (temperature)
  void eval_temp_mean(
    const double, //!< Height of the node from terrain
    double&       //!< Temperature source term to be populated
    );

  //! Inactive selector representing union of all the parts
  inline stk::mesh::Selector& inactive_selector() { return inactiveSelector_; }

private:
  ABLPostProcessingAlgorithm();
  ABLPostProcessingAlgorithm(const ABLPostProcessingAlgorithm&);

  //! Helper method that determines the part corresponding to a desired
  //! vertical level and ensures that part exists in the mesh database.
  void determine_part_names(
    std::vector<double>&, std::vector<std::string>&, bool, std::string&);

  //! Register velocity and temperature fields on the appropriate parts based
  //! on user input.
  void register_fields();

  //! Helper method to compute the statistics on z-planes
  void calc_stats();

  //! Reference to Realm
  Realm& realm_;

  //! Pointer to SpatialAveragingAlgorithm
  SpatialAveragingAlgorithm * spatialAvg_;
  
  //! Heights where velocity information is provided
  std::vector<double> heights_; // Array of shape [num_Uheights]

  //! Planar average velocity calculated on the surface [num_UHeights, 3]
  Array2D<double> UmeanCalc_;

  //! Planar average temperature calculated on the surface [num_THeights]
  std::vector<double> TmeanCalc_;

  //! stk::Transfer search methods
  std::string searchMethod_;
  //! stk::Transfer search tolerance
  double searchTolerance_;
  //! stk::Transfer search expansion factor
  double searchExpansionFactor_;

  //! Domains where velocity/temperature are averaged
  std::vector<std::string> fromTargetNames_;

  //! Part names
  std::vector<std::string> partNames_;
  std::unordered_set<std::string> allPartNames_;

  stk::mesh::PartVector allParts_;
  stk::mesh::Selector inactiveSelector_;

  Transfers* transfers_;

  //! Flag indicating whether to generate part names list for velocity field
  bool genPartList_;

  //! Format string specifier for generating velocity parts list
  std::string partFmt_;

  //! Write frequency for source term output
  int outputFreq_;

  //! Format string specifier indicating the file name for output. The
  //! specification takes one `%s` specifier that is used to populate Ux, Uy,
  //! Uz, T. Default is "abl_stats_%s.dat"
  std::string outFileFmt_;
};
}
}

#endif /* ABLPOSTPROCESSINGALGORITHM_H */
