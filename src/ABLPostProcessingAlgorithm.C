
#include "ABLPostProcessingAlgorithm.h"
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

ABLPostProcessingAlgorithm::ABLPostProcessingAlgorithm(Realm& realm, const YAML::Node& node)
  : realm_(realm),
    spatialAvg_(NULL),
    heights_(0),
    UmeanCalc_(0),
    TmeanCalc_(0),
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
    outFileFmt_("abl_%s_stats.dat")
{
  load(node);
}

ABLPostProcessingAlgorithm::ABLPostProcessingAlgorithm(Realm& realm, const YAML::Node& node, SpatialAveragingAlgorithm& spatialAvg)
  : realm_(realm),
    spatialAvg_(&spatialAvg),
    heights_(0),
    UmeanCalc_(0),
    TmeanCalc_(0),
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
    outFileFmt_("abl_%s_stats.dat")
{
    load(node);
}
    

ABLPostProcessingAlgorithm::~ABLPostProcessingAlgorithm()
{
  if (NULL != transfers_)
    delete transfers_;
}

void
ABLPostProcessingAlgorithm::load(const YAML::Node& node)
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
          "ABL Postprocessing: Mismatch between sizes of heights and target parts ");
      genPartList_ = false;
  } else if (node["target_part_format"]) { // Generate part names from printf
                                           // style string
      get_required<std::string>(node, "target_part_format", partFmt_);
      genPartList_ = true;
      partNames_.resize(nHeights);
  } else {
      throw std::runtime_error(
          "ABL Postprocessing: No target part(s) provided.");
  }

  if(spatialAvg_ == NULL) {
      spatialAvg_ = new SpatialAveragingAlgorithm(realm_, fromTargetNames_);
  }
  
  const int ndim = realm_.spatialDimension_;
  UmeanCalc_.resize(nHeights);
  for (size_t i = 0; i < nHeights; i++) {
      UmeanCalc_[i].resize(ndim);
  }
  TmeanCalc_.resize(nHeights);
  
}

void
ABLPostProcessingAlgorithm::setup()
{
  // Momentum sources
  determine_part_names(
    heights_, partNames_, genPartList_, partFmt_);
  // Register fields
  register_fields();
}

void
ABLPostProcessingAlgorithm::determine_part_names(
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
        "ABLPostProcessingAlgorithm::setup: Cannot find part " + partName);
    } else {
      // Add this part to the visited part list
      allPartNames_.insert(partName);
    }
  }
}

void
ABLPostProcessingAlgorithm::register_fields()
{
  stk::mesh::MetaData& meta = realm_.meta_data();
  int nDim = meta.spatial_dimension();
  int nStates = realm_.number_of_states();
  
  for (auto key : partNames_) {
    stk::mesh::Part* part = meta.get_part(key);
    VectorFieldType* vel = meta.get_field<VectorFieldType>
        (stk::topology::NODE_RANK, "velocity");
    spatialAvg_->register_part_field<VectorFieldType>(part, vel);
    ScalarFieldType* temp = meta.get_field<ScalarFieldType>(
        stk::topology::NODE_RANK, "temperature");
    spatialAvg_->register_part_field<ScalarFieldType>(part, temp);
  }
}

void
ABLPostProcessingAlgorithm::initialize()
{
  stk::mesh::MetaData& meta = realm_.meta_data();

  // We expect all parts to exist, so no creation step here

  // Add all parts to inactive selection
  for (auto key : allPartNames_) {
    stk::mesh::Part* part = meta.get_part(key);
    allParts_.push_back(part);
  }
  inactiveSelector_ = stk::mesh::selectUnion(allParts_);

  NaluEnv::self().naluOutputP0()
      << "ABL postprocessing active \n"
      << "\t Number of planes: " << heights_.size() << std::endl;

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
ABLPostProcessingAlgorithm::execute()
{
  calc_stats();
}

void
ABLPostProcessingAlgorithm::calc_stats()
{
  stk::mesh::MetaData& meta = realm_.meta_data();
  stk::mesh::BulkData& bulk = realm_.bulk_data();
  const int nDim = meta.spatial_dimension();
  const double dt = realm_.get_time_step();
  const double currTime = realm_.get_current_time();

  VectorFieldType* velocity =
    meta.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  ScalarFieldType* temperature =
    meta.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "temperature");

  const size_t nVarStats = 9 ;// <u'^2>, <v'^2>, <w'^2>, <u'v'>, <u'w'>, <v'w'>, <w'^3>, <theta'^2>, <w'theta'>
  const size_t numPlanes = heights_.size();
  // Sum(vectors) and number of nodes on this processor over all planes
  std::vector<double> sumVarStats(numPlanes * nVarStats, 0.0);
  std::vector<unsigned> numNodes(numPlanes, 0);
  // Global sum and nodes for computing global average
  std::vector<double> sumVarStatsGlobal(numPlanes * nVarStats, 0.0);
  std::vector<unsigned> totalNodes(numPlanes, 0);
  for (size_t ih = 0; ih < numPlanes; ih++) {
    stk::mesh::Part* part = meta.get_part(partNames_[ih]);
    spatialAvg_->eval_mean(velocity, part, UmeanCalc_[ih].data());
    spatialAvg_->eval_mean(temperature, part, &TmeanCalc_[ih]);
    
    const int ioff = ih * nVarStats;
    stk::mesh::Selector s_local_part(*part);
    const stk::mesh::BucketVector& node_buckets =
      bulk.get_buckets(stk::topology::NODE_RANK, s_local_part);

    // Calculate sum(velocity) for all nodes on this processor
    for (size_t ib = 0; ib < node_buckets.size(); ib++) {
      stk::mesh::Bucket& bukt = *node_buckets[ib];
      double* velField = stk::mesh::field_data(*velocity, bukt);
      double* tempField = stk::mesh::field_data(*temperature, bukt);

      for (size_t in = 0; in < bukt.size(); in++) {
        const int offset = in * nDim;
        sumVarStats[ioff + 0] += (velField[offset + 0] - UmeanCalc_[ih][0])*(velField[offset + 0] - UmeanCalc_[ih][0]); //<u'^2>
        sumVarStats[ioff + 1] += (velField[offset + 1] - UmeanCalc_[ih][1])*(velField[offset + 1] - UmeanCalc_[ih][1]); //<v'^2>
        sumVarStats[ioff + 2] += (velField[offset + 2] - UmeanCalc_[ih][2])*(velField[offset + 2] - UmeanCalc_[ih][2]); //<w'^2>
        sumVarStats[ioff + 3] += (velField[offset + 0] - UmeanCalc_[ih][0])*(velField[offset + 1] - UmeanCalc_[ih][1]); //<u'v'>
        sumVarStats[ioff + 4] += (velField[offset + 0] - UmeanCalc_[ih][0])*(velField[offset + 2] - UmeanCalc_[ih][2]); //<u'w'>
        sumVarStats[ioff + 5] += (velField[offset + 1] - UmeanCalc_[ih][1])*(velField[offset + 2] - UmeanCalc_[ih][2]); //<v'w'>
        sumVarStats[ioff + 6] += (velField[offset + 2] - UmeanCalc_[ih][2])*(velField[offset + 2] - UmeanCalc_[ih][2])*(velField[offset + 2] - UmeanCalc_[ih][2]); //<w'^3>
        sumVarStats[ioff + 7] += (tempField[in] - TmeanCalc_[ih])*(tempField[in] - TmeanCalc_[ih]); //<theta'^2>
        sumVarStats[ioff + 8] += (velField[offset + 2] - UmeanCalc_[ih][2])*(tempField[in] - TmeanCalc_[ih]); //<w'theta'>        
      }
      numNodes[ih] += bukt.size();
    }
  }

  // Assemble global sum and node count
  stk::all_reduce_sum(
    NaluEnv::self().parallel_comm(), sumVarStats.data(), sumVarStatsGlobal.data(),
    numPlanes * nDim);
  // Revisit this for area or volume weighted averaging.
  stk::all_reduce_sum(
    NaluEnv::self().parallel_comm(), numNodes.data(), totalNodes.data(),
    numPlanes);

  // // Compute spatial averages
  // for (size_t ih = 0; ih < numPlanes; ih++) {
  //   const size_t ioff = ih * nDim;
  //   for (int i = 0; i < nDim; i++) {
  //     UmeanCalc_[ih][i] = sumVarStatsGlobal[ioff + i] / totalNodes[ih];
  //   }
  // }

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
ABLPostProcessingAlgorithm::eval_vel_mean(
  const double zp, std::vector<double>& velMean)
{
  const int nDim = realm_.spatialDimension_;
  if (heights_.size() == 1) {
    // Constant source term throughout the domain
    for (int i = 0; i < nDim; i++) {
      velMean[i] = UmeanCalc_[i][0];
    }
  } else {
    // Linearly interpolate source term within the planes, maintain constant
    // source term above and below the heights provided
    for (int i = 0; i < nDim; i++) {
      utils::linear_interp(heights_, UmeanCalc_[i], zp, velMean[i]);
    }
  }
}

void
ABLPostProcessingAlgorithm::eval_temp_mean(const double zp, double& tempMean)
{
  if (heights_.size() == 1) {
    tempMean = TmeanCalc_[0];
  } else {
    utils::linear_interp(heights_, TmeanCalc_, zp, tempMean);
  }
}

} // namespace nalu
} // namespace sierra
