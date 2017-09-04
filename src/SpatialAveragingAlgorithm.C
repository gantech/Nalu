
#include "SpatialAveragingAlgorithm.h"
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

#include <boost/format.hpp>
#include <fstream>
#include <iostream>
#include <iomanip>

namespace sierra {
namespace nalu {

SpatialAveragingAlgorithm::SpatialAveragingAlgorithm(Realm& realm, const YAML::Node& node)
  : realm_(realm),
    heights_(0),
    searchMethod_("stk_kdtree"),
    searchTolerance_(1.0e-4),
    searchExpansionFactor_(1.5),
    fromTargetNames_(0),
    partNames_(),
    inactiveSelector_(),
    transfers_(0),
    genPartList_(false),
    nVectorAvg_(0),
    nScalarAvg_(0),
    vectorAvg_(NULL),
    scalarAvg_(NULL),
    partFmt_(""),
    outputFreq_(10),
    outFile_("spatial_averaging.dat")
{
  load(node);
}

SpatialAveragingAlgorithm::SpatialAveragingAlgorithm(Realm& realm, const std::vector<std::string>& fromTargetNames)
  : realm_(realm),
    heights_(0),
    searchMethod_("stk_kdtree"),
    searchTolerance_(1.0e-4),
    searchExpansionFactor_(1.5),
    fromTargetNames_(fromTargetNames),
    partNames_(),
    inactiveSelector_(),
    transfers_(0),
    genPartList_(false),
    nVectorAvg_(0),
    nScalarAvg_(0),
    vectorAvg_(NULL),
    scalarAvg_(NULL),
    partFmt_(""),
    outputFreq_(10),
    outFile_("")
{
}
    
SpatialAveragingAlgorithm::~SpatialAveragingAlgorithm()
{
  if (NULL != transfers_)
    delete transfers_;
}

void
SpatialAveragingAlgorithm::load(const YAML::Node& node)
{
  get_if_present(node, "search_method", searchMethod_, searchMethod_);
  get_if_present(node, "search_tolerance", searchTolerance_, searchTolerance_);
  get_if_present(node, "search_expansion_factor", searchExpansionFactor_,
                 searchExpansionFactor_);
  
  get_if_present(node, "output_frequency", outputFreq_, outputFreq_);
  get_if_present(node, "output_file", outFile_, outFile_);

  if (node["from_target_part"]) {
    const YAML::Node& fromParts = node["from_target_part"];
    if (fromParts.IsSequence()) {
      fromTargetNames_ = fromParts.as<std::vector<std::string>>();
    } else if (fromParts.IsScalar()) {
      fromTargetNames_.resize(1);
      fromTargetNames_[0] = fromParts.as<std::string>();
    }
  } else {
    throw std::runtime_error(
      "No from_target_part specified for ABL postprocessing function");
  }

  get_required<std::vector<double>>(node, "heights", heights_);
  auto nHeights = heights_.size();

  if (node["target_parts"]) { // Explicit target parts list provided
      get_required<std::vector<std::string>>(node, "target_parts", partNames_);
      ThrowAssertMsg(
          (nHeights == partNames_.size()),
          "Spatial Averaging: Mismatch between sizes of heights and target parts ");
      genPartList_ = false;
  } else if (node["target_part_format"]) { // Generate part names from printf
                                           // style string
      get_required<std::string>(node, "target_part_format", partFmt_);
      genPartList_ = true;
      partNames_.resize(nHeights);
  } else {
      throw std::runtime_error(
          "Spatial Averaging: No target part(s) provided.");
  }

  // extract the output variables
  const YAML::Node y_outputs = expect_sequence(node, "output_variables", false);
  if (y_outputs) {
      for (size_t ioutput = 0; ioutput < y_outputs.size(); ++ioutput) {
          const YAML::Node y_output = y_outputs[ioutput];
          // find the name, size and type
          const YAML::Node fieldNameNode = y_output["field_name"];
          const YAML::Node fieldSizeNode = y_output["field_size"];
    
          if ( !fieldNameNode ) 
              throw std::runtime_error("SpatialAveraging::load() Sorry, field name must be provided");
            
          if ( !fieldSizeNode ) 
              throw std::runtime_error("SpatialAveraging::load() Sorry, field size must be provided");
          
          fieldName_.push_back(fieldNameNode.as<std::string>());
          fieldSize_.push_back(fieldSizeNode.as<int>()) ;
      }
  }

  const int ndim = realm_.spatialDimension_;
  
}

void
SpatialAveragingAlgorithm::setup()
{
  // Momentum sources
  determine_part_names(
    heights_, partNames_, genPartList_, partFmt_);
  // Register fields
  register_fields();
}

void
SpatialAveragingAlgorithm::determine_part_names(
  std::vector<double>& heights,
  std::vector<std::string>& nameSet,
  bool flag,
  std::string& nameFmt)
{
  stk::mesh::MetaData& meta = realm_.meta_data();

  if (flag)
    for (size_t i = 0; i < heights.size(); i++)
      nameSet[i] = (boost::format(nameFmt) % heights[i]).str();

  for (auto partName : nameSet) {
    auto it = allPartNames_.find(partName);
    if (it != allPartNames_.end())
      continue;

    stk::mesh::Part* part = meta.get_part(partName);
    if (NULL == part) {
      // TODO: Need "nodeset" creation capability. Integrate with
      // DataProbePostProcessing to minimize code duplication.
      throw std::runtime_error(
        "SpatialAveragingAlgorithm::setup: Cannot find part " + partName);
    } else {
      // Add this part to the visited part list
      allPartNames_.insert(partName);
    }
  }
}

void
SpatialAveragingAlgorithm::register_fields()
{
  stk::mesh::MetaData& meta = realm_.meta_data();
  int nDim = meta.spatial_dimension();
  int nStates = realm_.number_of_states();

  for (auto key : allPartNames_) {
    stk::mesh::Part* part = meta.get_part(key);
    VectorFieldType& coords = meta.declare_field<VectorFieldType>(
      stk::topology::NODE_RANK, "coordinates");
    stk::mesh::put_field(coords, *part, nDim);
  }

  for (std::size_t iPart=0; iPart < partNames_.size(); iPart++) {
      stk::mesh::Part* part = meta.get_part(partNames_[iPart]);
      for(std::size_t jField=0; jField < fieldName_.size(); jField++) {
          if (fieldSize_[jField] == 1) {
              ScalarFieldType& scalarField = meta.declare_field<ScalarFieldType>(
                  stk::topology::NODE_RANK, fieldName_[jField], nStates);
              stk::mesh::put_field(scalarField, *part, nDim);
              scalarAvgMap_.insert( { {partNames_[iPart], fieldName_[jField]}, nScalarAvg_ } );
              scalarIO_.push_back(true);
	      nScalarAvg_++;
          } else if (fieldSize_[jField] == 3) {
              VectorFieldType& vecField = meta.declare_field<VectorFieldType>(
                  stk::topology::NODE_RANK, fieldName_[jField], nStates);
              stk::mesh::put_field(vecField, *part, nDim);
              vectorAvgMap_.insert( { {partNames_[iPart], fieldName_[jField]}, nVectorAvg_ } );
              vectorIO_.push_back(true);
	      nVectorAvg_++;
          }
      }
  }

}

void
SpatialAveragingAlgorithm::initialize()
{
  stk::mesh::MetaData& meta = realm_.meta_data();
  int nDim = meta.spatial_dimension();

  vectorAvg_.resize(nVectorAvg_);
  for(size_t i=0; i < nVectorAvg_; i++)
    vectorAvg_[i].resize(nDim);
  scalarAvg_.resize(nScalarAvg_);

  // Add all parts to inactive selection
  for (auto key : allPartNames_) {
    stk::mesh::Part* part = meta.get_part(key);
    allParts_.push_back(part);
  }
  inactiveSelector_ = stk::mesh::selectUnion(allParts_);
  create_transfers();

  NaluEnv::self().naluOutputP0()
      << "Spatial Averaging active \n"
      << "\t Number of planes: " << allPartNames_.size() << std::endl;

  // Prepare output files to dump sources when computed during precursor phase
  if (( NaluEnv::self().parallel_rank() == 0 )) {
    std::fstream spAvgFile;
    spAvgFile.open(outFile_.c_str(), std::fstream::out);
    spAvgFile.close();
  }
}

void
SpatialAveragingAlgorithm::create_transfers()
{
  transfers_ = new Transfers(*realm_.root());

  for(auto key: vectorAvgMap_) {
      const std::pair<std::string, std::string> partFieldPair = key.first;
      populate_transfer_data(partFieldPair.second, partFieldPair.first);
  }
  
  for(auto key: scalarAvgMap_) {
      const std::pair<std::string, std::string> partFieldPair = key.first;
      populate_transfer_data(partFieldPair.second, partFieldPair.first);
  }

  transfers_->initialize();
}

void
SpatialAveragingAlgorithm::populate_transfer_data(
    const std::string & fieldName, const std::string & partName)
{
  stk::mesh::MetaData& meta = realm_.meta_data();

  Transfer* theXfer = new Transfer(*transfers_);
  theXfer->name_ = "spatialAveraging_xfer_" + partName + "_" + fieldName;
  theXfer->fromRealm_ = &realm_;
  theXfer->toRealm_ = &realm_;
  theXfer->searchMethodName_ = searchMethod_;
  theXfer->searchTolerance_ = searchTolerance_;
  theXfer->searchExpansionFactor_ = searchExpansionFactor_;

  for (auto key : fromTargetNames_) {
    stk::mesh::Part* part = meta.get_part(key);
    theXfer->fromPartVec_.push_back(part);
  }
  
  stk::mesh::Part* part = meta.get_part(partName);
  theXfer->toPartVec_.push_back(part);

  theXfer->transferVariablesPairName_.push_back(
    std::make_pair(fieldName, fieldName));
  transfers_->transferVector_.push_back(theXfer);
}

void
SpatialAveragingAlgorithm::execute()
{
  // Map fields from fluidRealm onto averaging planes
  transfers_->execute();
  calc_scalar_averages();
  calc_vector_averages();
}

void
SpatialAveragingAlgorithm::calc_vector_averages()
{
  stk::mesh::MetaData& meta = realm_.meta_data();
  stk::mesh::BulkData& bulk = realm_.bulk_data();
  const int nDim = meta.spatial_dimension();
  const double currTime = realm_.get_current_time();

  // Sum(vectors) and number of nodes on this processor over all planes
  std::vector<double> sumVect(nVectorAvg_ * nDim, 0.0);
  std::vector<unsigned> numNodes(nVectorAvg_, 0);
  // Global sum and nodes for computing global average
  std::vector<double> sumVectGlobal(nVectorAvg_ * nDim, 0.0);
  std::vector<unsigned> totalNodes(nVectorAvg_, 0);

  for(auto key: vectorAvgMap_) {
      const std::pair<std::string, std::string> partFieldPair = key.first;
      const size_t ivAvg = key.second;
      const int ioff = ivAvg * nDim;
      VectorFieldType* vectField =
          meta.get_field<VectorFieldType>(stk::topology::NODE_RANK, partFieldPair.second);
      stk::mesh::Part* part = meta.get_part(partFieldPair.first);
      stk::mesh::Selector s_local_part(*part);
      const stk::mesh::BucketVector& node_buckets =
          bulk.get_buckets(stk::topology::NODE_RANK, s_local_part);
      // Calculate sum(vectors) for all nodes on this processor
      for (size_t ib = 0; ib < node_buckets.size(); ib++) {
          stk::mesh::Bucket& bukt = *node_buckets[ib];
          double* vectFieldP = stk::mesh::field_data(*vectField, bukt);
          
          for (size_t in = 0; in < bukt.size(); in++) {
              const int offset = in * nDim;
              for (int i = 0; i < nDim; i++) {
                  sumVect[ioff + i] += vectFieldP[offset + i];
	      }
          }
          numNodes[ivAvg] += bukt.size();
      }
  }

  // Assemble global sum and node count
  stk::all_reduce_sum(
    NaluEnv::self().parallel_comm(), sumVect.data(), sumVectGlobal.data(),
    nVectorAvg_ * nDim);
  // Revisit this for area or volume weighted averaging.
  stk::all_reduce_sum(
    NaluEnv::self().parallel_comm(), numNodes.data(), totalNodes.data(),
    nVectorAvg_);

  // Compute spatial averages 
  for(auto key: vectorAvgMap_) {
      const size_t ivAvg = key.second;
      for (int i = 0; i < nDim; i++) 
	vectorAvg_[ivAvg][i] = sumVectGlobal[ivAvg*nDim + i] / totalNodes[ivAvg];

  }

  // Write averages to file if necessary
  const int tcount = realm_.get_time_step_count();
  if (( NaluEnv::self().parallel_rank() == 0 ) &&
      ( tcount % outputFreq_ == 0)) {
      std::fstream spAvgFile;
      spAvgFile.open(outFile_.c_str(), std::fstream::app);
      spAvgFile << "Time: " << std::setw(12) << currTime << std::endl;
      for(auto key: vectorAvgMap_) {
          const size_t ivAvg = key.second;
          if(vectorIO_[ivAvg]) {
              spAvgFile << "Part: " << (key.first).first << " Field: " << (key.first).second ;
              for (int i = 0; i < nDim; i++)                  
                  spAvgFile << std::setprecision(6)
                            << std::setw(15)
                            << vectorAvg_[ivAvg][i];
              spAvgFile << std::endl ;
          }
      }
      spAvgFile << std::endl;
      spAvgFile.close();
  }
 
}

void
SpatialAveragingAlgorithm::calc_scalar_averages()
{
  stk::mesh::MetaData& meta = realm_.meta_data();
  stk::mesh::BulkData& bulk = realm_.bulk_data();
  const int nDim = meta.spatial_dimension();
  const double currTime = realm_.get_current_time();

  // Sum(vectors) and number of nodes on this processor over all planes
  std::vector<double> sumScal(nScalarAvg_, 0.0);
  std::vector<unsigned> numNodes(nScalarAvg_, 0);
  // Global sum and nodes for computing global average
  std::vector<double> sumScalGlobal(nScalarAvg_, 0.0);
  std::vector<unsigned> totalNodes(nScalarAvg_, 0);

  for(auto key: scalarAvgMap_) {
      const std::pair<std::string, std::string> partFieldPair = key.first;
      const size_t isAvg = key.second;
      ScalarFieldType* scalField =
          meta.get_field<ScalarFieldType>(stk::topology::NODE_RANK, partFieldPair.second);
      stk::mesh::Part* part = meta.get_part(partFieldPair.first);
      stk::mesh::Selector s_local_part(*part);
      const stk::mesh::BucketVector& node_buckets =
          bulk.get_buckets(stk::topology::NODE_RANK, s_local_part);
      // Calculate sum(scalars) for all nodes on this processor
      for (size_t ib = 0; ib < node_buckets.size(); ib++) {
          stk::mesh::Bucket& bukt = *node_buckets[ib];
          double* scalFieldP = stk::mesh::field_data(*scalField, bukt);
          
          for (size_t in = 0; in < bukt.size(); in++) 
              sumScal[isAvg] += scalFieldP[in];
          numNodes[isAvg] += bukt.size();
      }
  }

  // Assemble global sum and node count
  stk::all_reduce_sum(
    NaluEnv::self().parallel_comm(), sumScal.data(), sumScalGlobal.data(),
    nScalarAvg_);
  // Revisit this for area or volume weighted averaging.
  stk::all_reduce_sum(
    NaluEnv::self().parallel_comm(), numNodes.data(), totalNodes.data(),
    nScalarAvg_);

  // Compute spatial averages
  for(auto key: scalarAvgMap_) {
      const size_t isAvg = key.second;
      scalarAvg_[isAvg] = sumScalGlobal[isAvg] / totalNodes[isAvg];
  }

  // Write averages to file if necessary
  const int tcount = realm_.get_time_step_count();
  if (( NaluEnv::self().parallel_rank() == 0 ) &&
      ( tcount % outputFreq_ == 0)) {
      std::fstream spAvgFile;
      spAvgFile.open(outFile_.c_str(), std::fstream::app);
      spAvgFile << "Time: " << std::setw(12) << currTime << std::endl;
      for(auto key: scalarAvgMap_) {
          const size_t isAvg = key.second;
          if(vectorIO_[isAvg]) {
              spAvgFile << "Part: " << (key.first).first << " Field: " << (key.first).second ;
              spAvgFile << std::setprecision(6)
                        << std::setw(15)
                        << scalarAvg_[isAvg]
                        << std::endl ;
      }
      spAvgFile << std::endl;
      spAvgFile.close();
  }
}

}


} // namespace nalu
} // namespace sierra
