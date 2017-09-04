
#ifndef SPATIALAVERAGINGALGORITHM_H
#define SPATIALAVERAGINGALGORITHM_H

#include "NaluParsing.h"
#include "FieldTypeDef.h"
#include "Realm.h"
#include "xfer/Transfer.h"
#include "xfer/Transfers.h"
#include "utils/LinearInterpolation.h"

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

#include <stk_io/IossBridge.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

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

  //! Store average value in maps from pair(partName, field_id) to id of average value array
  size_t nVectorAvg_;
  size_t nScalarAvg_;
  std::vector<std::vector<double>> vectorAvg_;
  std::vector<double> scalarAvg_;
  std::map<std::pair<std::string, std::string>, size_t > vectorAvgMap_;
  std::map<std::pair<std::string, std::string>, size_t > scalarAvgMap_;
  
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

template<typename FieldType>
void SpatialAveragingAlgorithm::register_part_field(
    stk::mesh::Part* part, FieldType* field)
{
    stk::mesh::MetaData& meta = realm_.meta_data();
    int nDim = meta.spatial_dimension();
    int nStates = realm_.number_of_states();

    auto it = allPartNames_.find(part->name());
    if (it != allPartNames_.end()) {
        // Part already exists. Nothing to do here
    } else {
        allPartNames_.insert(part->name());
        VectorFieldType& coords = meta.declare_field<VectorFieldType>(
            stk::topology::NODE_RANK, "coordinates");
        stk::mesh::put_field(coords, *part, nDim);
    }
    
    std::pair<std::string, std::string> partFieldPair(part->name(), field->name());
    if(field->field_array_rank()) {
        if(vectorAvgMap_.find(partFieldPair) != vectorAvgMap_.end()) {
            // Part - Field combination already exists. Nothing to do here  
        } else {
            VectorFieldType& vecField = meta.declare_field<VectorFieldType>(
                stk::topology::NODE_RANK, field->name(), nStates);
            stk::mesh::put_field(vecField, *part, nDim);
            vectorAvgMap_.insert( { {part->name(), field->name()}, nVectorAvg_ } );
	    nVectorAvg_++;
        }
    } else {
        if(scalarAvgMap_.find(partFieldPair) != scalarAvgMap_.end()) {
            // Part - Field combination already exists. Nothing to do here  
        } else {
            ScalarFieldType& scalarField = meta.declare_field<ScalarFieldType>(
                stk::topology::NODE_RANK, field->name());
            stk::mesh::put_field(scalarField, *part);
            scalarAvgMap_.insert( { {part->name(), field->name()}, nScalarAvg_ } );
	    nScalarAvg_++;
        }
    }
    
}

template<typename FieldType>
void SpatialAveragingAlgorithm::eval_mean(
    FieldType* field, stk::mesh::Part* part, double * avgValue)
{
    const int nDim = realm_.spatialDimension_;

    std::pair<std::string, std::string> partFieldPair(part->name(), field->name());
    if(field->field_array_rank()) {

      auto ivMap = vectorAvgMap_.find(partFieldPair);
      if( ivMap != vectorAvgMap_.end()) {
        const size_t ivAvg = ivMap->second;
        for(int i=0; i<nDim; i++)
	  avgValue[i] = vectorAvg_[ivAvg][i];

      } else {
	NaluEnv::self().naluOutput() << "Part field pair not found. Part:" << part->name() << ", Field: " << field->name() << std::endl ;
	NaluEnv::self().naluOutput() << "Available pairs are :" << std::endl ;
	for (auto vavgPairs: vectorAvgMap_) {
	  NaluEnv::self().naluOutputP0() << (vavgPairs.first).first << " " << (vavgPairs.first).second <<  std::endl ;
	}
      }
    } else {
      
      auto isMap = scalarAvgMap_.find(partFieldPair);
      if( isMap != scalarAvgMap_.end()) {
        const size_t isAvg = isMap->second;
        *avgValue = scalarAvg_[isAvg];
      } else {
	NaluEnv::self().naluOutput() << "Part field pair not found. Part:" << part->name() << ", Field: " << field->name() << std::endl ;
	NaluEnv::self().naluOutput() << "Available pairs are :" << std::endl ;
	for (auto savgPairs: scalarAvgMap_) {
	  NaluEnv::self().naluOutput() << (savgPairs.first).first << " " << (savgPairs.first).second <<  std::endl ;
	}
      }
    }
    
}


}
}

#endif /* SPATIALAVERAGINGALGORITHM_H */
