/**
 * @file descartes_robot_sampler.hpp
 * @brief Tesseract Descartes Kinematics Sampler Implementation
 *
 * @author Levi Armstrong
 * @date April 18, 2018
 * @version TODO
 * @bug No known bugs
 *
 * @copyright Copyright (c) 2017, Southwest Research Institute
 *
 * @par License
 * Software License Agreement (Apache License)
 * @par
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0
 * @par
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#ifndef TESSERACT_MOTION_PLANNERS_DESCARTES_ROBOT_SAMPLER_HPP
#define TESSERACT_MOTION_PLANNERS_DESCARTES_ROBOT_SAMPLER_HPP

#include <tesseract_common/macros.h>
TESSERACT_COMMON_IGNORE_WARNINGS_PUSH
#include <console_bridge/console.h>
#include <Eigen/Geometry>
#include <vector>
TESSERACT_COMMON_IGNORE_WARNINGS_POP

#include <tesseract_motion_planners/descartes/descartes_robot_sampler.h>
#include <tesseract_kinematics/core/utils.h>

namespace tesseract_planning
{
template <typename FloatType>
DescartesRobotSampler<FloatType>::DescartesRobotSampler(std::string target_working_frame,
                                                        const Eigen::Isometry3d& target_pose,
                                                        PoseSamplerFn target_pose_sampler,
                                                        tesseract_kinematics::KinematicGroup::ConstPtr manip,
                                                        DescartesCollision::Ptr collision,
                                                        std::string tcp_frame,
                                                        const Eigen::Isometry3d& tcp_offset,
                                                        bool allow_collision,
                                                        DescartesVertexEvaluator::Ptr is_valid,
                                                        bool use_redundant_joint_solutions)
  : target_working_frame_(std::move(target_working_frame))
  , target_pose_(target_pose)
  , target_pose_sampler_(std::move(target_pose_sampler))
  , manip_(std::move(manip))
  , collision_(std::move(collision))
  , tcp_frame_(std::move(tcp_frame))
  , tcp_offset_(tcp_offset)
  , allow_collision_(allow_collision)
  , dof_(static_cast<int>(manip_->numJoints()))
  , ik_seed_(Eigen::VectorXd::Zero(dof_))
  , is_valid_(std::move(is_valid))
  , use_redundant_joint_solutions_(use_redundant_joint_solutions)
{
  if (!allow_collision_ && !collision_)
    throw std::runtime_error("Collision checker must not be a nullptr if collisions are not allowed during planning");
}

template <typename FloatType>
std::vector<descartes_light::StateSample<FloatType>> DescartesRobotSampler<FloatType>::sample() const
{
  // Generate all possible Cartesian poses
  tesseract_common::VectorIsometry3d target_poses = target_pose_sampler_(target_pose_);

  bool found_ik_sol = false;
  std::stringstream error_string_stream;

  auto console_bridge_loglevel = console_bridge::getLogLevel();

  // Generate the IK solutions for those poses
  std::vector<descartes_light::StateSample<FloatType>> samples;
  //  for (const auto& pose : target_poses)
  for (std::size_t i = 0; i < target_poses.size(); i++)
  {
    const auto& pose = target_poses[i];

    // Get the transformation to the kinematic tip link
    Eigen::Isometry3d target_pose = pose * tcp_offset_.inverse();

    // Solve IK (TODO Should tcp_offset be stored in KinGroupIKInput?)
    tesseract_kinematics::KinGroupIKInput ik_input(target_pose, target_working_frame_, tcp_frame_);
    tesseract_kinematics::IKSolutions ik_solutions = manip_->calcInvKin({ ik_input }, ik_seed_);

    if (ik_solutions.empty())
      continue;

    tesseract_collision::ContactTrajectoryResults traj_contacts(manip_->getJointNames(), ik_solutions.size());

    found_ik_sol = true;

    // Check each individual joint solution
    for (std::size_t j = 0; j < ik_solutions.size(); j++)
    {
      const auto& sol = ik_solutions[j];

      if ((is_valid_ != nullptr) && !(*is_valid_)(sol))
        continue;

      auto state = std::make_shared<descartes_light::State<FloatType>>(sol.cast<FloatType>());
      if (allow_collision_ && collision_ == nullptr)
      {
        samples.push_back(descartes_light::StateSample<FloatType>{ state, static_cast<FloatType>(0.0) });
      }
      else if (!allow_collision_)
      {
        if (console_bridge_loglevel == console_bridge::LogLevel::CONSOLE_BRIDGE_LOG_DEBUG)
        {
          tesseract_collision::ContactResultMap coll_results = collision_->detailed_validate(sol);
          if (!coll_results.empty())
          {
            tesseract_collision::ContactTrajectoryStepResults step_contacts(j, sol, sol, 1);
            tesseract_collision::ContactTrajectorySubstepResults substep_contacts(1, sol);
            substep_contacts.contacts = coll_results;
            step_contacts.substeps[0] = substep_contacts;
            traj_contacts.steps[j] = step_contacts;
          }
          else
            samples.push_back(descartes_light::StateSample<FloatType>{ state, 0.0 });
        }
        else if (collision_->validate(sol))
          samples.push_back(descartes_light::StateSample<FloatType>{ state, 0.0 });
      }
      else
      {
        const FloatType cost = static_cast<FloatType>(collision_->distance(sol));
        samples.push_back(descartes_light::StateSample<FloatType>{ state, cost });
      }
    }

    if (console_bridge_loglevel == console_bridge::LogLevel::CONSOLE_BRIDGE_LOG_DEBUG)
    {
      error_string_stream << "For sample " << i << " the target position is:" << std::endl
                          << target_pose.matrix() << std::endl;
      error_string_stream << ik_solutions.size()
                          << " IK solutions were found, with a collision summary of:" << std::endl;
      error_string_stream << traj_contacts.trajectoryCollisionResultsTable().str();
    }
  }

  if (samples.empty())
  {
    Eigen::Isometry3d first_pose = target_poses.front();
    Eigen::Isometry3d first_pose_offset = first_pose * tcp_offset_.inverse();
    std::stringstream ss;
    ss << "Descartes vertex failure: ";
    if (!found_ik_sol)
      ss << "No IK solutions were found. ";
    else
      ss << "All IK solutions found were in collision or invalid. ";
    ss << target_poses.size() << " samples tried from target pose" << std::endl;
    ss << "\ttarget pose translation: (" << target_pose_.translation().x() << ", " << target_pose_.translation().y()
       << ", " << target_pose_.translation().z() << ")" << std::endl;
    ;
    ss << "\tworking frame: '" << target_working_frame_ << "', tcp: '" << tcp_frame_ << "'" << std::endl;
    if (found_ik_sol)
      ss << error_string_stream.str();
    if (console_bridge_loglevel == console_bridge::LogLevel::CONSOLE_BRIDGE_LOG_DEBUG)
      std::cout << ss.str();
    return samples;
  }

  if (allow_collision_)
  {
    // Sort state samples in descending order by distance from nearest collision (i.e. state cost)
    std::sort(samples.begin(), samples.end(), [](const auto& a, const auto& b) { return a.cost > b.cost; });

    // Convert the distance into a cost and normalize
    if (samples.size() > 1)
    {
      const FloatType max_dist = samples.front().cost;
      const FloatType min_dist = samples.back().cost;
      const FloatType range = max_dist - min_dist;
      if (range > std::numeric_limits<FloatType>::epsilon())
      {
        std::for_each(samples.begin(), samples.end(), [&min_dist, &range](auto& sample) {
          sample.cost = static_cast<FloatType>(1.0) - (sample.cost - min_dist) / range;
        });
      }
      else
      {
        std::for_each(samples.begin(), samples.end(), [](auto& sample) { sample.cost = 0.0; });
      }
    }
  }

  // Generate the redundant solutions
  if (use_redundant_joint_solutions_)
  {
    const Eigen::MatrixX2d& limits = manip_->getLimits().joint_limits;
    std::vector<Eigen::Index> redundancy_capable_joints = manip_->getRedundancyCapableJointIndices();
    std::vector<descartes_light::StateSample<FloatType>> redundant_samples;
    for (const descartes_light::StateSample<FloatType>& sample : samples)
    {
      const auto redundant_solutions =
          tesseract_kinematics::getRedundantSolutions<FloatType>(*(sample.state), limits, redundancy_capable_joints);

      // Add the redundant samples with the same cost as the nominal sample
      std::transform(redundant_solutions.begin(),
                     redundant_solutions.end(),
                     std::back_inserter(redundant_samples),
                     [&sample](const auto& sol) {
                       auto state = std::make_shared<descartes_light::State<FloatType>>(sol.template cast<FloatType>());
                       return descartes_light::StateSample<FloatType>{ state, sample.cost };
                     });
    }

    // Combine the nominal samples with the redundant samples
    samples.insert(samples.end(), redundant_samples.begin(), redundant_samples.end());
  }

  return samples;
}

template <typename FloatType>
std::ostream& DescartesRobotSampler<FloatType>::format(std::ostream& out) const
{
  auto tr = this->target_pose_.translation();
  auto rot = Eigen::Quaternion<double>(this->target_pose_.rotation());
  out << "{ [x=" << tr.x() << ", y=" << tr.y() << ", z=" << tr.z() << "], [x=" << rot.x() << ", y=" << rot.y()
      << ", z=" << rot.z() << ", w=" << rot.w() << "] }";
  return out;
}

}  // namespace tesseract_planning

#endif  // TESSERACT_MOTION_PLANNERS_DESCARTES_ROBOT_SAMPLER_HPP
