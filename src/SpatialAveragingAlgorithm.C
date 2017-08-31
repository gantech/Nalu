
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
    partFmt_(""),
    outputFreq_(10),
    outFileFmt_("spatial_averaging_%s.dat")
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
    partFmt_(""),
    outputFreq_(10),
    outFileFmt_("")
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
  get_if_present(node, "output_format", outFileFmt_, outFileFmt_);

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
              scalarAvg_.insert( { {partNames_[iPart], fieldName_[jField]}, double(0.0) } );
          } else if (fieldSize_[jField] == 3) {
              VectorFieldType& vecField = meta.declare_field<VectorFieldType>(
                  stk::topology::NODE_RANK, fieldName_[jField], nStates);
              stk::mesh::put_field(vecField, *part, nDim);
              std::vector<double> zeroVector;
              zeroVector.resize(nDim);
              for(int i=0; i<nDim; i++) zeroVector[i] = 0.0;
              vectorAvg_.insert( { {partNames_[iPart], fieldName_[jField]}, zeroVector } );
          }
      }
  }

}

template<typename FieldType>
void SpatialAveragingAlgorithm::register_part_field(stk::mesh::Part* part, FieldType* field) {

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
        if(vectorAvg_.find(partFieldPair)) {
            // Part - Field combination already exists. Nothing to do here  
        } else {
            VectorFieldType& vecField = meta.declare_field<VectorFieldType>(
                stk::topology::NODE_RANK, field->name(), nStates);
            stk::mesh::put_field(vecField, *part, nDim);
            std::vector<double> zeroVector;
            zeroVector.resize(nDim);
            for(int i=0; i<nDim; i++) zeroVector[i] = 0.0;
            vectorAvg_.insert( { {part->name(), field->name()}, zeroVector } );            
        }
    } else {
        if(scalarAvg_.find(partFieldPair)) {
            // Part - Field combination already exists. Nothing to do here  
        } else {
            ScalarFieldType& scalarField = meta.declare_field<ScalarFieldType>(
                stk::topology::NODE_RANK, field->name(), nStates);
            stk::mesh::put_field(scalarField, *part, nDim);
            scalarAvg_.insert( { {part->name(), field->name()}, 0.0 } );
        }
    }
    
}

void
SpatialAveragingAlgorithm::initialize()
{
  stk::mesh::MetaData& meta = realm_.meta_data();

  // We expect all parts to exist, so no creation step here

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
    std::string uxname((boost::format(outFileFmt_)%"Ux").str());
    std::string uyname((boost::format(outFileFmt_)%"Uy").str());
    std::string uzname((boost::format(outFileFmt_)%"Uz").str());
    std::fstream uxFile, uyFile, uzFile;
    uxFile.open(uxname.c_str(), std::fstream::out);
    uyFile.open(uyname.c_str(), std::fstream::out);
    uzFile.open(uzname.c_str(), std::fstream::out);

    uxFile << "# Time, " ;
    uyFile << "# Time, " ;
    uzFile << "# Time, " ;
    for (size_t ih = 0; ih < heights_.size(); ih++) {
      uxFile << heights_[ih] << ", ";
      uyFile << heights_[ih] << ", ";
      uzFile << heights_[ih] << ", ";
    }
    uxFile << std::endl ;
    uyFile << std::endl ;
    uzFile << std::endl ;
    uxFile.close();
    uyFile.close();
    uzFile.close();
  }
}

void
SpatialAveragingAlgorithm::create_transfers()
{
  transfers_ = new Transfers(*realm_.root());

  for(auto key: vectorAvg_) {
      const std::pair<std::string, std::string> partFieldPair = key.first;
      populate_transfer_data(partFieldPair.second, partFieldPair.first);
  }
  
  for(auto key: scalarAvg_) {
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
  const double dt = realm_.get_time_step();
  const double currTime = realm_.get_current_time();

  const size_t numParts = vectorAvg_.size();
  // Sum(vectors) and number of nodes on this processor over all planes
  std::vector<double> sumVect(numParts * nDim, 0.0);
  std::vector<unsigned> numNodes(numParts, 0);
  // Global sum and nodes for computing global average
  std::vector<double> sumVectGlobal(numParts * nDim, 0.0);
  std::vector<unsigned> totalNodes(numParts, 0);

  size_t iPart=0;
  for(auto key: vectorAvg_) {
      const std::pair<std::string, std::string> partFieldPair = key.first;
      const int ioff = iPart * nDim;
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
              for (int i = 0; i < nDim; i++)
                  sumVect[ioff + i] += vectFieldP[offset + i];
          }
          numNodes[iPart] += bukt.size();
      }
      iPart++;
  }

  // Assemble global sum and node count
  stk::all_reduce_sum(
    NaluEnv::self().parallel_comm(), sumVect.data(), sumVectGlobal.data(),
    numParts * nDim);
  // Revisit this for area or volume weighted averaging.
  stk::all_reduce_sum(
    NaluEnv::self().parallel_comm(), numNodes.data(), totalNodes.data(),
    numParts);

  // Compute spatial averages
  iPart=0;
  for(auto key: vectorAvg_) {
      std::vector<double> vAvg;
      vAvg.resize(nDim);
      for (int i = 0; i < nDim; i++)
          vAvg[i] = sumVectGlobal[iPart*nDim + i] / totalNodes[iPart];
      key.second = vAvg;
      iPart++;
  }

  // const int tcount = realm_.get_time_step_count();
  // if (( NaluEnv::self().parallel_rank() == 0 ) &&
  //     ( tcount % outputFreq_ == 0)) {
  //   std::string uxname((boost::format(outFileFmt_)%"Ux").str());
  //   std::string uyname((boost::format(outFileFmt_)%"Uy").str());
  //   std::string uzname((boost::format(outFileFmt_)%"Uz").str());
  //   std::fstream uxFile, uyFile, uzFile;
  //   uxFile.open(uxname.c_str(), std::fstream::app);
  //   uyFile.open(uyname.c_str(), std::fstream::app);
  //   uzFile.open(uzname.c_str(), std::fstream::app);

  //   uxFile << std::setw(12) << currTime << ", ";
  //   uyFile << std::setw(12) << currTime << ", ";
  //   uzFile << std::setw(12) << currTime << ", ";
  //   for (size_t ih = 0; ih < heights_.size(); ih++) {
  //     uxFile << std::setprecision(6)
  //            << std::setw(15)
  //            << UmeanCalc_[0][ih] << ", ";
  //     uyFile << std::setprecision(6)
  //            << std::setw(15)
  //            << UmeanCalc_[1][ih] << ", ";
  //     uzFile << std::setprecision(6)
  //            << std::setw(15)
  //            << UmeanCalc_[2][ih] << ", ";
  //   }
  //   uxFile << std::endl;
  //   uyFile << std::endl;
  //   uzFile << std::endl;
  //   uxFile.close();
  //   uyFile.close();
  //   uzFile.close();
  // }
  
}

void
SpatialAveragingAlgorithm::calc_scalar_averages()
{
  stk::mesh::MetaData& meta = realm_.meta_data();
  stk::mesh::BulkData& bulk = realm_.bulk_data();
  const int nDim = meta.spatial_dimension();
  const double dt = realm_.get_time_step();
  const double currTime = realm_.get_current_time();

  const size_t numParts = scalarAvg_.size();
  // Sum(vectors) and number of nodes on this processor over all planes
  std::vector<double> sumScal(numParts, 0.0);
  std::vector<unsigned> numNodes(numParts, 0);
  // Global sum and nodes for computing global average
  std::vector<double> sumScalGlobal(numParts, 0.0);
  std::vector<unsigned> totalNodes(numParts, 0);

  size_t iPart=0;
  for(auto key: scalarAvg_) {
      const std::pair<std::string, std::string> partFieldPair = key.first;
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
              sumScal[iPart] += scalFieldP[in];
          numNodes[iPart] += bukt.size();
      }
      iPart++;
  }

  // Assemble global sum and node count
  stk::all_reduce_sum(
    NaluEnv::self().parallel_comm(), sumScal.data(), sumScalGlobal.data(),
    numParts);
  // Revisit this for area or volume weighted averaging.
  stk::all_reduce_sum(
    NaluEnv::self().parallel_comm(), numNodes.data(), totalNodes.data(),
    numParts);

  // Compute spatial averages
  iPart=0;
  for(auto key: scalarAvg_) {
      key.second = sumScalGlobal[iPart] / totalNodes[iPart];
      iPart++;
  }

  // const int tcount = realm_.get_time_step_count();
  // if (( NaluEnv::self().parallel_rank() == 0 ) &&
  //     ( tcount % outputFreq_ == 0)) {
  //   std::string uxname((boost::format(outFileFmt_)%"Ux").str());
  //   std::string uyname((boost::format(outFileFmt_)%"Uy").str());
  //   std::string uzname((boost::format(outFileFmt_)%"Uz").str());
  //   std::fstream uxFile, uyFile, uzFile;
  //   uxFile.open(uxname.c_str(), std::fstream::app);
  //   uyFile.open(uyname.c_str(), std::fstream::app);
  //   uzFile.open(uzname.c_str(), std::fstream::app);

  //   uxFile << std::setw(12) << currTime << ", ";
  //   uyFile << std::setw(12) << currTime << ", ";
  //   uzFile << std::setw(12) << currTime << ", ";
  //   for (size_t ih = 0; ih < heights_.size(); ih++) {
  //     uxFile << std::setprecision(6)
  //            << std::setw(15)
  //            << UmeanCalc_[0][ih] << ", ";
  //     uyFile << std::setprecision(6)
  //            << std::setw(15)
  //            << UmeanCalc_[1][ih] << ", ";
  //     uzFile << std::setprecision(6)
  //            << std::setw(15)
  //            << UmeanCalc_[2][ih] << ", ";
  //   }
  //   uxFile << std::endl;
  //   uyFile << std::endl;
  //   uzFile << std::endl;
  //   uxFile.close();
  //   uyFile.close();
  //   uzFile.close();
  // }
  
}

template<typename FieldType>
void SpatialAveragingAlgorithm::eval_mean(
    FieldType* field, stk::mesh::Part* part, double * avgValue)
{
    const int nDim = realm_.spatialDimension_;

    std::pair<std::string, std::string> partFieldPair(part->name(), field->name());
    if(field->field_array_rank()) {
        std::vector<double> mapValue = vectorAvg_.find(partFieldPair)->second;
        for(int i=0; i<nDim; i++)
            avgValue[i] = mapValue[i];
    } else {
        std::vector<double> mapValue = scalarAvg_.find(partFieldPair)->second;
        *avgValue = mapValue[0];
    }
    
}

} // namespace nalu
} // namespace sierra
