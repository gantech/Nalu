#include "ABLPostProcessingAlgorithm.h"
#include "SpatialAveragingAlgorithm.h"
#include "Realm.h"
#include <master_element/MasterElement.h>
#include "xfer/Transfer.h"
#include "xfer/Transfers.h"
#include "utils/LinearInterpolation.h"
#include "Plane2D.h"

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
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
    indepSpatAvg_(true),
    heights_(0),
    SFSstressMeanCalc_(0),
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

ABLPostProcessingAlgorithm::ABLPostProcessingAlgorithm(Realm& realm, const YAML::Node& node, SpatialAveragingAlgorithm* spatialAvg)
  : realm_(realm),
    spatialAvg_(spatialAvg),
    indepSpatAvg_(false),
    heights_(0),
    SFSstressMeanCalc_(0),
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
  if(indepSpatAvg_ && (NULL != spatialAvg_))
    delete spatialAvg_;

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

  if (node["abl_wall_parts"]) {
      const YAML::Node& ablWallParts = node["abl_wall_parts"];
      if (ablWallParts.IsSequence()) {
          ablWallNames_ = ablWallParts.as<std::vector<std::string>>();
      } else if (ablWallParts.IsScalar()) {
          ablWallNames_.resize(1);
          ablWallNames_[0] = ablWallParts.as<std::string>();
      }
  } else {
      throw std::runtime_error(
          "No abl_wall_parts specified for ABL postprocessing function");
  }

  get_if_present(node, "generate_parts", generateParts_, generateParts_);
  if (generateParts_) {
      quadVertices_ = node["quad_vertices"].as<std::vector<std::vector<double>>>();
      ThrowAssertMsg(quadVertices_.size() == 4,
                     "ABLPostProcessing::load: Invalid data encountered when parsing 'quad_vertices'");
      for (auto crd: quadVertices_) {
          ThrowAssertMsg(crd.size() == 2,
                         "ABLPostProcessing::load: Invalid data encountered when parsing 'quad_vertices'");
      }
      nx_ = node["nx"].as<int>();
      ny_ = node["ny"].as<int>();
  }
  
  if(spatialAvg_ == NULL) {
      spatialAvg_ = new SpatialAveragingAlgorithm(realm_, fromTargetNames_);
  }

  const int ndim = realm_.spatialDimension_;
  SFSstressMeanCalc_.resize(nHeights);
  UmeanCalc_.resize(nHeights);
  varCalc_.resize(nHeights);
  for (size_t i = 0; i < nHeights; i++) {
      SFSstressMeanCalc_[i].resize(6);
      UmeanCalc_[i].resize(ndim);
      varCalc_[i].resize(nVarStats_);
  }
  TmeanCalc_.resize(nHeights);

}

void
ABLPostProcessingAlgorithm::setup()
{

  if(indepSpatAvg_)
    spatialAvg_->setup() ;
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
  stk::mesh::BulkData& bulk = realm_.bulk_data();

  if (flag)
    for (size_t i = 0; i < heights.size(); i++)
      nameSet[i] = (boost::format(nameFmt) % heights[i]).str();

  for (size_t i=0; i < heights.size(); i++) {
    auto partName = nameSet[i];
    auto it = allPartNames_.find(partName);
    if (it != allPartNames_.end())
      continue;

    stk::mesh::Part* part = meta.get_part(partName);
    if (NULL == part) {
        if (!generateParts_)
            throw std::runtime_error(
                "ABLPostProcessingAlgorithm::setup: Cannot find part " + partName);
        else {
            planeGenerators_.emplace_back(
                new Plane2D(meta, bulk, partName));
            Plane2D& plane = *planeGenerators_[planeGenerators_.size() - 1];
            auto& vertices = plane.vertices();

            for (int j=0; j<4; j++) {
                vertices[j][0] = quadVertices_[j][0];
                vertices[j][1] = quadVertices_[j][1];
                vertices[j][2] = heights[i];
            }
            plane.set_dimensions(nx_, ny_);
            plane.setup();
        }
    } 
    // Add this part to the visited part list
    allPartNames_.insert(partName);
  }

  for (auto partName : ablWallNames_) {
      stk::mesh::Part* part = meta.get_part(partName);
      if (NULL == part) {
          throw std::runtime_error(
              "ABLPostProcessingAlgorithm::setup: Cannot find part " + partName);
      } else {

          ablWallPartVec_.push_back(part);
      }
  }

}

void
ABLPostProcessingAlgorithm::register_fields()
{
  stk::mesh::MetaData& meta = realm_.meta_data();

  for (auto key : partNames_) {
    stk::mesh::Part* part = meta.get_part(key);
    GenericFieldType* sfsStress = meta.get_field<GenericFieldType>
        (stk::topology::NODE_RANK, "sfs_stress");
    spatialAvg_->register_part_field<GenericFieldType>(part, sfsStress, 6);
    VectorFieldType* vel = meta.get_field<VectorFieldType>
        (stk::topology::NODE_RANK, "velocity");
    spatialAvg_->register_part_field<VectorFieldType>(part, vel, 3);
    ScalarFieldType* temp = meta.get_field<ScalarFieldType>(
        stk::topology::NODE_RANK, "temperature");
    spatialAvg_->register_part_field<ScalarFieldType>(part, temp, 1);
  }
}

void
ABLPostProcessingAlgorithm::initialize()
{
  stk::mesh::MetaData& meta = realm_.meta_data();

  for (auto& plane: planeGenerators_)
    plane->initialize();

  if(indepSpatAvg_)
    spatialAvg_->initialize() ;

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
    std::string uMeanName((boost::format(outFileFmt_)%"U").str());
    std::string tMeanName((boost::format(outFileFmt_)%"T").str());
    std::string varName((boost::format(outFileFmt_)%"Var").str());
    std::string uTauName((boost::format(outFileFmt_)%"uTau").str());
    std::fstream uMeanFile, tMeanFile, varFile, uTauFile;
    uMeanFile.open(uMeanName.c_str(), std::fstream::out);
    tMeanFile.open(tMeanName.c_str(), std::fstream::out);
    varFile.open(varName.c_str(), std::fstream::out);
    uTauFile.open(uTauName.c_str(), std::fstream::out);
    uMeanFile << "# Time, " ;
    tMeanFile << "# Time, " ;
    varFile << "# Time, " ;
    uTauFile << "# Time, " ;
    for (size_t ih = 0; ih < heights_.size(); ih++) {
      for (size_t iDim = 0; iDim < 3; iDim++)
          uMeanFile << heights_[ih] << ", ";
      tMeanFile << heights_[ih] << ", ";
      for (size_t iVarStat = 0; iVarStat < nVarStats_; iVarStat++)
          varFile << heights_[ih] << ", ";
      for (size_t iSFSstress = 0; iSFSstress < 6; iSFSstress++)
          varFile << heights_[ih] << ", ";
    }
    uMeanFile << std::endl ;
    tMeanFile << std::endl ;
    varFile << std::endl ;
    uTauFile << std::endl ;
    uMeanFile.close();
    tMeanFile.close();
    varFile.close() ;
    uTauFile.close();
  }
}

void
ABLPostProcessingAlgorithm::execute()
{
  if(indepSpatAvg_)
    spatialAvg_->execute() ;
  calc_stats();
  calc_utau();
}

void
ABLPostProcessingAlgorithm::calc_stats()
{
  stk::mesh::MetaData& meta = realm_.meta_data();
  stk::mesh::BulkData& bulk = realm_.bulk_data();
  const double currTime = realm_.get_current_time();

  GenericFieldType* sfs_stress =
      meta.get_field<GenericFieldType>(stk::topology::NODE_RANK, "sfs_stress");
  VectorFieldType* velocity =
    meta.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  ScalarFieldType* temperature =
    meta.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "temperature");

  const size_t numPlanes = heights_.size();
  // Sum(vectors) and number of nodes on this processor over all planes
  std::vector<double> sumVarStats(numPlanes * nVarStats_, 0.0);
  std::vector<unsigned> numNodes(numPlanes, 0);
  // Global sum and nodes for computing global average
  std::vector<double> sumVarStatsGlobal(numPlanes * nVarStats_, 0.0);
  std::vector<unsigned> totalNodes(numPlanes, 0);
  for (size_t ih = 0; ih < numPlanes; ih++) {
    stk::mesh::Part* part = meta.get_part(partNames_[ih]);

    spatialAvg_->eval_mean(sfs_stress, part, &(SFSstressMeanCalc_[ih][0]), 6);
    spatialAvg_->eval_mean(velocity, part, &(UmeanCalc_[ih][0]), 3);
    spatialAvg_->eval_mean(temperature, part, &TmeanCalc_[ih], 1);

    const int ioff = ih * nVarStats_;
    stk::mesh::Selector s_local_part(*part);
    const stk::mesh::BucketVector& node_buckets =
      bulk.get_buckets(stk::topology::NODE_RANK, s_local_part);

    // Calculate sum(velocity) for all nodes on this processor
    for (size_t ib = 0; ib < node_buckets.size(); ib++) {
      stk::mesh::Bucket& bukt = *node_buckets[ib];
      double* velField = stk::mesh::field_data(*velocity, bukt);
      double* tempField = stk::mesh::field_data(*temperature, bukt);

      for (size_t in = 0; in < bukt.size(); in++) {
        const int offset = in * 3;
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
    numPlanes * nVarStats_);
  // Revisit this for area or volume weighted averaging.
  stk::all_reduce_sum(
    NaluEnv::self().parallel_comm(), numNodes.data(), totalNodes.data(),
    numPlanes);

  // Compute variances and covariances
  for (size_t ih = 0; ih < numPlanes; ih++) {
    const size_t ioff = ih * nVarStats_;
    for (size_t i = 0; i < nVarStats_; i++) {
      varCalc_[ih][i] = sumVarStatsGlobal[ioff + i] / totalNodes[ih];
    }
  }

  const int tcount = realm_.get_time_step_count();
  if (( NaluEnv::self().parallel_rank() == 0 ) &&
      ( tcount % outputFreq_ == 0)) {
    std::string uMeanName((boost::format(outFileFmt_)%"U").str());
    std::string tMeanName((boost::format(outFileFmt_)%"T").str());
    std::string varName((boost::format(outFileFmt_)%"Var").str());
    std::fstream uMeanFile, tMeanFile, varFile;
    uMeanFile.open(uMeanName.c_str(), std::fstream::app);
    tMeanFile.open(tMeanName.c_str(), std::fstream::app);
    varFile.open(varName.c_str(), std::fstream::app);
    uMeanFile << std::setw(12) << currTime << ", ";
    tMeanFile << std::setw(12) << currTime << ", ";
    varFile << std::setw(12) << currTime << ", ";
    for (size_t ih = 0; ih < heights_.size(); ih++) {
      for (size_t iDim = 0; iDim < 3; iDim++)
          uMeanFile << std::setprecision(6)
                    << std::setw(15)
                    << UmeanCalc_[ih][iDim] << ", ";
      tMeanFile << std::setprecision(6)
                << std::setw(15)
                << TmeanCalc_[ih] << ", ";
      for (size_t iVarStat = 0; iVarStat < nVarStats_; iVarStat++)
          varFile << std::setprecision(6)
                  << std::setw(15)
                  << varCalc_[ih][iVarStat] << ", ";
      for (size_t iSFSstress = 0; iSFSstress < 6; iSFSstress++)
          varFile << std::setprecision(6)
                  << std::setw(15)
                  << SFSstressMeanCalc_[ih][iSFSstress] << ", ";

    }
    uMeanFile << std::endl;
    tMeanFile << std::endl;
    varFile << std::endl ;
    uMeanFile.close();
    tMeanFile.close();
    varFile.close() ;
  }

}


void
ABLPostProcessingAlgorithm::calc_utau()
{

    std::array<double, 2> uTauAreaSumLocal{0.0, 0.0}; // First value hold uTau * area, second value holds area. To get average utau, sum both quantities over all processes and divide one by another.
    std::array<double, 2> uTauAreaSumGlobal{0.0, 0.0};

    stk::mesh::MetaData & meta_data = realm_.meta_data();
    const double currTime = realm_.get_current_time();

    GenericFieldType* exposedAreaVec = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "exposed_area_vector");
    GenericFieldType* wallFrictionVelocityBip = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "wall_friction_velocity_bip");
    stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
        &stk::mesh::selectUnion(ablWallPartVec_);
    stk::mesh::BucketVector const& face_buckets =
        realm_.get_buckets( meta_data.side_rank(), s_locally_owned_union );
    for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
          ib != face_buckets.end() ; ++ib ) {
        stk::mesh::Bucket & b = **ib ;

        // face master element
        MasterElement *meFC = sierra::nalu::get_surface_master_element(b.topology());
        const int numScsBip = meFC->numIntPoints_;

        const stk::mesh::Bucket::size_type length   = b.size();

        for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

            // get face
            stk::mesh::Entity face = b[k];

            // pointer to face data
            const double * areaVec = stk::mesh::field_data(*exposedAreaVec, face);
            double *wallFVBip = stk::mesh::field_data(*wallFrictionVelocityBip, face);

            // loop over face nodes
            for ( int ip = 0; ip < numScsBip; ++ip ) {

                const int offSetAveraVec = ip*3;
                double aMag = 0.0;
                for ( int j = 0; j < 3; ++j ) {
                    const double axj = areaVec[offSetAveraVec+j];
                    aMag += axj*axj;
                }
                aMag = std::sqrt(aMag);

                uTauAreaSumLocal[0] += wallFVBip[ip] * aMag ;
                uTauAreaSumLocal[1] += aMag ;
            }
        }
    }

    stk::all_reduce_sum(
        NaluEnv::self().parallel_comm(), uTauAreaSumLocal.data(), uTauAreaSumGlobal.data(), 2);

    utauCalc_ = uTauAreaSumGlobal[0]/uTauAreaSumGlobal[1] ;

    const int tcount = realm_.get_time_step_count();
    if (( NaluEnv::self().parallel_rank() == 0 ) &&
        ( tcount % outputFreq_ == 0)) {
        std::string uTauName((boost::format(outFileFmt_)%"uTau").str());
        std::fstream uTauFile;
        uTauFile.open(uTauName.c_str(), std::fstream::app);
        uTauFile << std::setw(12) << currTime << ", ";
        uTauFile << std::setprecision(6)
                 << std::setw(15)
                 << utauCalc_ ;
        uTauFile << std::endl;
        uTauFile.close();
    }

}

void
ABLPostProcessingAlgorithm::eval_vel_mean(
  const double zp, std::vector<double>& velMean)
{
  if (heights_.size() == 1) {
    // Constant source term throughout the domain
    for (int i = 0; i < 3; i++) {
      velMean[i] = UmeanCalc_[i][0];
    }
  } else {
    // Linearly interpolate source term within the planes, maintain constant
    // source term above and below the heights provided
    for (int i = 0; i < 3; i++) {
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
