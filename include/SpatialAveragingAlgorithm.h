
#ifndef SPATIALAVERAGINGALGORITHM_H
#define SPATIALAVERAGINGALGORITHM_H

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

/**
 * \brief Spatial averaging
 *
 * This class parses the user inputs and provides an spatial average within Nalu.
 * The SpatialAveraging capability is turned on by the presence of a sub-section
 * titled `spatial_averaging` within the Realm section of the Nalu input file.
 *
 * ```
 *   spatial_averaging:
 *     search_method: stk_kdtree
 *     search_tolerance: 0.0001
 *     search_expansion_factor: 1.5
 *     output_frequency: 10
 *     from_target_part: [Unspecified-2-HEX]
 *     target_part: [Unspecified-2-QUAD, Unspecified-3-QUAD]
 *     target_part_format: "zplane_%.1f"
 *     heights: [80.0, 160.0]
 *     variables:
 *        - field_name: velocity
 *          field_size: 3
 *
 *        - field_name: temperature
 *          field_size: 1
 * ```
 *
 */
class SpatialAveragingAlgorithm
{
public:
  template <typename T>
  using Array2D = std::vector<std::vector<T>>;

  SpatialAveragingAlgorithm(Realm&, const YAML::Node&);
  SpatialAveragingAlgorithm(Realm&, const std::vector<std::string>& fromTargetNames);
  
  ~SpatialAveragingAlgorithm();

  //! Parse input file for user options and initialize
  void load(const YAML::Node&);

  //! Setup ABL postprocessing (steps before mesh creation)
  void setup();

  //! Register part-field combination to average
  template<typename FieldType>
      void register_part_field(stk::mesh::Part* part, FieldType* field);
  
  //! Initialize ABL postprocessing (steps after mesh creation)
  void initialize();

  //! Execute field transfers, compute planar averaging of fields
  //! on parts
  void execute();

  //! Evaluate the spatial average at
  template<typename FieldType>
      void eval_mean(FieldType * field, stk::mesh::Part * part, double * avgValue) ;
      
  //! Inactive selector representing union of all the parts
  inline stk::mesh::Selector& inactive_selector() { return inactiveSelector_; }

private:
  SpatialAveragingAlgorithm();
  SpatialAveragingAlgorithm(const SpatialAveragingAlgorithm&);

  //! Helper method that determines the part corresponding to a desired
  //! vertical level and ensures that part exists in the mesh database.
  void determine_part_names(
    std::vector<double>&, std::vector<std::string>&, bool, std::string&);

  //! Register fields on the appropriate parts based on user input.
  void register_fields();

  //! Create transfer that handles mapping of fields from
  //! fluidRealm to the planar nodesets.
  void create_transfers();

  void populate_transfer_data(const std::string &, const std::string &);

  //! Calculate averages
  void calc_scalar_averages();
  void calc_vector_averages();
  
  //! Reference to Realm
  Realm& realm_;
  
  //! Heights where velocity information is provided
  std::vector<double> heights_; // Array of shape [num_heights]

  //! Store average value in maps from pair(partName, field_id) to average value
  std::map<std::pair<std::string, std::string>, std::vector<double> > vectorAvg_;
  std::map<std::pair<std::string, std::string>, double > scalarAvg_;
  
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

  std::vector<std::string> fieldName_;
  std::vector<int> fieldSize_;
  
  Transfers* transfers_;

  //! Flag indicating whether to generate part names list for velocity field
  bool genPartList_;

  //! Format string specifier for generating velocity parts list
  std::string partFmt_;

  //! width for output
  int w_;

 //! Write frequency for source term output
  int outputFreq_;

  //! Format string specifier indicating the file name for output. The
  //! specification takes one `%s` specifier. Default is
  //! "spatial_average_%s.dat"
  std::string outFileFmt_;
};
}
}

#endif /* SPATIALAVERAGINGALGORITHM_H */
