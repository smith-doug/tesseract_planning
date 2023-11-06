/**
 * @file ompl_motion_planner.cpp
 * @brief Tesseract OMPL motion planner
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

#include <tesseract_common/macros.h>
TESSERACT_COMMON_IGNORE_WARNINGS_PUSH
#include <console_bridge/console.h>
#include <ompl/base/goals/GoalState.h>
#include <ompl/base/goals/GoalStates.h>
#include <ompl/tools/multiplan/ParallelPlan.h>
#include <ompl/geometric/PathGeometric.h>
TESSERACT_COMMON_IGNORE_WARNINGS_POP

#include <tesseract_environment/utils.h>

#include <tesseract_motion_planners/planner_utils.h>

#include <tesseract_motion_planners/ompl/ompl_motion_planner_status_category.h>
#include <tesseract_motion_planners/ompl/ompl_motion_planner.h>
#include <tesseract_motion_planners/ompl/continuous_motion_validator.h>
#include <tesseract_motion_planners/ompl/discrete_motion_validator.h>
#include <tesseract_motion_planners/ompl/profile/ompl_default_plan_profile.h>
#include <tesseract_motion_planners/ompl/weighted_real_vector_state_sampler.h>
#include <tesseract_motion_planners/core/utils.h>

#include <tesseract_command_language/command_language.h>
#include <tesseract_command_language/utils/utils.h>

namespace tesseract_planning
{
bool checkStartState(const ompl::base::ProblemDefinitionPtr& prob_def,
                     const Eigen::Ref<const Eigen::VectorXd>& state,
                     const OMPLStateExtractor& extractor)
{
  if (!(prob_def->getStartStateCount() >= 1))
    return false;

  for (unsigned i = 0; i < prob_def->getStartStateCount(); ++i)
    if (extractor(prob_def->getStartState(i)).isApprox(state, 1e-5))
      return true;

  return false;
}

bool checkGoalState(const ompl::base::ProblemDefinitionPtr& prob_def,
                    const Eigen::Ref<const Eigen::VectorXd>& state,
                    const OMPLStateExtractor& extractor)
{
  ompl::base::GoalPtr goal = prob_def->getGoal();
  if (goal->getType() == ompl::base::GoalType::GOAL_STATE)
    return extractor(prob_def->getGoal()->as<ompl::base::GoalState>()->getState()).isApprox(state, 1e-5);

  if (goal->getType() == ompl::base::GoalType::GOAL_STATES)
  {
    auto* goal_states = prob_def->getGoal()->as<ompl::base::GoalStates>();
    for (unsigned i = 0; i < goal_states->getStateCount(); ++i)
      if (extractor(goal_states->getState(i)).isApprox(state, 1e-5))
        return true;
  }
  else
  {
    CONSOLE_BRIDGE_logWarn("checkGoalStates: Unsupported Goal Type!");
    return true;
  }
  return false;
}

/** @brief Construct a basic planner */
OMPLMotionPlanner::OMPLMotionPlanner(std::string name)
  : name_(std::move(name)), status_category_(std::make_shared<const OMPLMotionPlannerStatusCategory>(name_))
{
  if (name_.empty())
    throw std::runtime_error("OMPLMotionPlanner name is empty!");
}

const std::string& OMPLMotionPlanner::getName() const { return name_; }

bool OMPLMotionPlanner::terminate()
{
  CONSOLE_BRIDGE_logWarn("Termination of ongoing optimization is not implemented yet");
  return false;
}

tesseract_common::StatusCode OMPLMotionPlanner::solve(const PlannerRequest& request,
                                                      PlannerResponse& response,
                                                      bool verbose) const
{
  if (!checkUserInput(request))  // NOLINT
  {
    response.status =
        tesseract_common::StatusCode(OMPLMotionPlannerStatusCategory::ErrorInvalidInput, status_category_);
    return response.status;
  }
  std::vector<OMPLProblem::Ptr> problem;
  if (request.data)
  {
    problem = *std::static_pointer_cast<std::vector<OMPLProblem::Ptr>>(request.data);
  }
  else
  {
    try
    {
      problem = createProblems(request);
    }
    catch (std::exception& e)
    {
      CONSOLE_BRIDGE_logError("OMPLPlanner failed to generate problem: %s.", e.what());
      response.status =
          tesseract_common::StatusCode(OMPLMotionPlannerStatusCategory::ErrorInvalidInput, status_category_);
      return response.status;
    }

    response.data = std::make_shared<std::vector<OMPLProblem::Ptr>>(problem);
  }

  // If the verbose set the log level to debug.
  if (verbose)
    console_bridge::setLogLevel(console_bridge::LogLevel::CONSOLE_BRIDGE_LOG_DEBUG);

  /// @todo: Need to expand this to support multiple motion plans leveraging taskflow
  for (auto& p : problem)
  {
    p->simple_setup->setup();
    auto parallel_plan = std::make_shared<ompl::tools::ParallelPlan>(p->simple_setup->getProblemDefinition());

    for (const auto& planner : p->planners)
      parallel_plan->addPlanner(planner->create(p->simple_setup->getSpaceInformation()));

    ompl::base::PlannerStatus status;
    if (!p->optimize)
    {
      // Solve problem. Results are stored in the response
      // Disabling hybridization because there is a bug which will return a trajectory that starts at the end state
      // and finishes at the end state.
      status = parallel_plan->solve(p->planning_time, 1, static_cast<unsigned>(p->max_solutions), false);
    }
    else
    {
      ompl::time::point end = ompl::time::now() + ompl::time::seconds(p->planning_time);
      const ompl::base::ProblemDefinitionPtr& pdef = p->simple_setup->getProblemDefinition();
      while (ompl::time::now() < end)
      {
        // Solve problem. Results are stored in the response
        // Disabling hybridization because there is a bug which will return a trajectory that starts at the end state
        // and finishes at the end state.
        ompl::base::PlannerStatus localResult =
            parallel_plan->solve(std::max(ompl::time::seconds(end - ompl::time::now()), 0.0),
                                 1,
                                 static_cast<unsigned>(p->max_solutions),
                                 false);
        if (localResult)
        {
          if (status != ompl::base::PlannerStatus::EXACT_SOLUTION)
            status = localResult;

          if (!pdef->hasOptimizationObjective())
          {
            CONSOLE_BRIDGE_logDebug("Terminating early since there is no optimization objective specified");
            break;
          }

          ompl::base::Cost obj_cost = pdef->getSolutionPath()->cost(pdef->getOptimizationObjective());
          CONSOLE_BRIDGE_logDebug("Motion Objective Cost: %f", obj_cost.value());

          if (pdef->getOptimizationObjective()->isSatisfied(obj_cost))
          {
            CONSOLE_BRIDGE_logDebug("Terminating early since solution path satisfies the optimization objective");
            break;
          }

          if (pdef->getSolutionCount() >= static_cast<std::size_t>(p->max_solutions))
          {
            CONSOLE_BRIDGE_logDebug("Terminating early since %u solutions were generated", p->max_solutions);
            break;
          }
        }
      }
    }

    if (status != ompl::base::PlannerStatus::EXACT_SOLUTION)
    {
      response.status = tesseract_common::StatusCode(OMPLMotionPlannerStatusCategory::ErrorFailedToFindValidSolution,
                                                     status_category_);
      return response.status;
    }

    if (p->simplify)
    {
      p->simple_setup->simplifySolution();
    }
    else
    {
      // Interpolate the path if it shouldn't be simplified and there are currently fewer states than requested
      auto num_output_states = static_cast<unsigned>(p->n_output_states);
      if (p->simple_setup->getSolutionPath().getStateCount() < num_output_states)
      {
        p->simple_setup->getSolutionPath().interpolate(num_output_states);
      }
      else
      {
        // Now try to simplify the trajectory to get it under the requested number of output states
        // The interpolate function only executes if the current number of states is less than the requested
        p->simple_setup->simplifySolution();
        if (p->simple_setup->getSolutionPath().getStateCount() < num_output_states)
          p->simple_setup->getSolutionPath().interpolate(num_output_states);
      }
    }
  }

  // Flatten the results to make them easier to process
  response.results = request.seed;
  std::vector<std::reference_wrapper<Instruction>> results_flattened =
      flattenProgramToPattern(response.results, request.instructions);
  std::vector<std::reference_wrapper<const Instruction>> instructions_flattened = flattenProgram(request.instructions);

  std::size_t instructions_idx = 0;  // Index for each input instruction

  // Handle the start instruction
  const auto& plan_instruction = instructions_flattened.at(0).get().as<PlanInstruction>();
  if (plan_instruction.isStart())
  {
    const auto& p = problem[0];

    // Get the results
    tesseract_common::TrajArray trajectory = p->getTrajectory();

    // Enforce limits
    for (Eigen::Index i = 0; i < trajectory.rows(); i++)
      tesseract_common::enforcePositionLimits(trajectory.row(i), p->manip->getLimits().joint_limits);

    assert(checkStartState(p->simple_setup->getProblemDefinition(), trajectory.row(0), p->extractor));
    assert(checkGoalState(p->simple_setup->getProblemDefinition(), trajectory.bottomRows(1).transpose(), p->extractor));

    // Copy the start instruction
    assert(instructions_idx == 0);
    assert(isMoveInstruction(results_flattened[0].get()));
    auto& move_instruction = results_flattened[0].get().as<MoveInstruction>();
    move_instruction.getWaypoint().as<StateWaypoint>().position = trajectory.row(0);
    instructions_idx++;
  }

  // Loop over remaining instructions
  std::size_t prob_idx = 0;
  for (; instructions_idx < instructions_flattened.size(); instructions_idx++)
  {
    if (isPlanInstruction(instructions_flattened.at(instructions_idx).get()))
    {
      const auto& p = problem[prob_idx];

      // Get the results
      tesseract_common::TrajArray trajectory = p->getTrajectory();

      assert(checkStartState(p->simple_setup->getProblemDefinition(), trajectory.row(0), p->extractor));
      assert(
          checkGoalState(p->simple_setup->getProblemDefinition(), trajectory.bottomRows(1).transpose(), p->extractor));

      // Loop over the flattened results and add them to response if the input was a plan instruction
      auto& move_instructions = results_flattened[instructions_idx].get().as<CompositeInstruction>();
      // Adjust result index to align final point since start instruction is already handled
      Eigen::Index result_index = trajectory.rows() - static_cast<Eigen::Index>(move_instructions.size());
      for (auto& instruction : move_instructions)
        instruction.as<MoveInstruction>().getWaypoint().as<StateWaypoint>().position = trajectory.row(result_index++);

      // Increment the problem
      prob_idx++;
    }
  }

  response.status = tesseract_common::StatusCode(OMPLMotionPlannerStatusCategory::SolutionFound, status_category_);
  return response.status;
}

void OMPLMotionPlanner::clear() { parallel_plan_ = nullptr; }

MotionPlanner::Ptr OMPLMotionPlanner::clone() const { return std::make_shared<OMPLMotionPlanner>(name_); }

bool OMPLMotionPlanner::checkUserInput(const PlannerRequest& request)
{
  // Check that parameters are valid
  if (request.env == nullptr)
  {
    CONSOLE_BRIDGE_logError("In TrajOptPlannerUniversalConfig: env is a required parameter and has not been set");
    return false;
  }

  if (request.instructions.empty())
  {
    CONSOLE_BRIDGE_logError("TrajOptPlannerUniversalConfig requires at least one instruction");
    return false;
  }

  return true;
}

OMPLProblem::Ptr createOMPLSubProblem(const PlannerRequest& request,
                                      const tesseract_kinematics::JointGroup::ConstPtr& manip)
{
  auto sub_prob = std::make_unique<OMPLProblem>();
  sub_prob->env = request.env;
  sub_prob->env_state = request.env_state;
  sub_prob->manip = manip;
  sub_prob->contact_checker = request.env->getDiscreteContactManager();
  sub_prob->contact_checker->setCollisionObjectsTransform(request.env_state.link_transforms);
  sub_prob->contact_checker->setActiveCollisionObjects(manip->getActiveLinkNames());
  return sub_prob;
}

std::vector<OMPLProblem::Ptr> OMPLMotionPlanner::createProblems(const PlannerRequest& request) const
{
  std::vector<OMPLProblem::Ptr> problem;
  tesseract_kinematics::JointGroup::Ptr manip;

  // Assume all the plan instructions have the same manipulator as the composite
  assert(!request.instructions.getManipulatorInfo().empty());

  const ManipulatorInfo& composite_mi = request.instructions.getManipulatorInfo();

  if (composite_mi.manipulator.empty())
    throw std::runtime_error("OMPL, manipulator is empty!");

  try
  {
    tesseract_kinematics::KinematicGroup::Ptr kin_group;
    std::string error_msg;
    if (composite_mi.manipulator_ik_solver.empty())
    {
      kin_group = request.env->getKinematicGroup(composite_mi.manipulator);
      error_msg = "Failed to find kinematic group for manipulator '" + composite_mi.manipulator + "'";
    }
    else
    {
      kin_group = request.env->getKinematicGroup(composite_mi.manipulator, composite_mi.manipulator_ik_solver);
      error_msg = "Failed to find kinematic group for manipulator '" + composite_mi.manipulator + "' with solver '" +
                  composite_mi.manipulator_ik_solver + "'";
    }

    if (kin_group == nullptr)
    {
      CONSOLE_BRIDGE_logError("%s", error_msg.c_str());
      throw std::runtime_error(error_msg);
    }

    manip = kin_group;
  }
  catch (...)
  {
    manip = request.env->getJointGroup(composite_mi.manipulator);
  }

  if (!manip)
    throw std::runtime_error("Failed to get joint/kinematic group: " + composite_mi.manipulator);

  std::vector<std::string> joint_names = manip->getJointNames();
  std::vector<std::string> active_link_names = manip->getActiveLinkNames();

  // Check and make sure it does not contain any composite instruction
  for (const auto& instruction : request.instructions)
    if (isCompositeInstruction(instruction))
      throw std::runtime_error("OMPL planner does not support child composite instructions.");

  int index = 0;
  Waypoint start_waypoint{ NullWaypoint() };
  Instruction placeholder_instruction{ NullInstruction() };
  const Instruction* start_instruction = nullptr;
  if (request.instructions.hasStartInstruction())
  {
    assert(isPlanInstruction(request.instructions.getStartInstruction()));
    start_instruction = &(request.instructions.getStartInstruction());
    if (isPlanInstruction(*start_instruction))
    {
      const auto& temp = start_instruction->as<PlanInstruction>();
      assert(temp.isStart());
      start_waypoint = temp.getWaypoint();
    }
    else
    {
      throw std::runtime_error("OMPL DefaultProblemGenerator: Unsupported start instruction type!");
    }
  }
  else
  {
    Eigen::VectorXd current_jv = request.env_state.getJointValues(joint_names);
    StateWaypoint swp(joint_names, current_jv);

    MoveInstruction temp_move(swp, MoveInstructionType::START);
    placeholder_instruction = temp_move;
    start_instruction = &placeholder_instruction;
    start_waypoint = swp;
  }

  // Transform plan instructions into ompl problem
  for (std::size_t i = 0; i < request.instructions.size(); ++i)
  {
    const auto& instruction = request.instructions[i];
    if (isPlanInstruction(instruction))
    {
      assert(isPlanInstruction(instruction));
      const auto& plan_instruction = instruction.as<PlanInstruction>();

      assert(isCompositeInstruction(request.seed[i]));
      const auto& seed_composite = request.seed[i].as<tesseract_planning::CompositeInstruction>();

      // Get Plan Profile
      std::string profile = plan_instruction.getProfile();
      profile = getProfileString(name_, profile, request.plan_profile_remapping);
      auto cur_plan_profile =
          getProfile<OMPLPlanProfile>(name_, profile, *request.profiles, std::make_shared<OMPLDefaultPlanProfile>());
      cur_plan_profile = applyProfileOverrides(name_, profile, cur_plan_profile, plan_instruction.profile_overrides);
      if (!cur_plan_profile)
        throw std::runtime_error("OMPLMotionPlannerDefaultConfig: Invalid profile");

      /** @todo Should check that the joint names match the order of the manipulator */
      OMPLProblem::Ptr sub_prob = createOMPLSubProblem(request, manip);
      cur_plan_profile->setup(*sub_prob);
      sub_prob->n_output_states = static_cast<int>(seed_composite.size()) + 1;

      if (plan_instruction.isLinear())
      {
        /** @todo Add support for linear motion to ompl planner */
        if (isCartesianWaypoint(plan_instruction.getWaypoint()) || isJointWaypoint(plan_instruction.getWaypoint()) ||
            isStateWaypoint(start_waypoint))
        {
          // TODO Currently skipping linear moves until SE3 motion planning is implemented.
          problem.push_back(nullptr);
          ++index;
        }
        else
        {
          throw std::runtime_error("OMPLMotionPlannerDefaultConfig: unknown waypoint type");
        }
      }
      else if (plan_instruction.isFreespace())
      {
        if (isJointWaypoint(plan_instruction.getWaypoint()) || isStateWaypoint(plan_instruction.getWaypoint()))
        {
          assert(checkJointPositionFormat(joint_names, plan_instruction.getWaypoint()));
          const Eigen::VectorXd& cur_position = getJointPosition(plan_instruction.getWaypoint());
          cur_plan_profile->applyGoalStates(
              *sub_prob, cur_position, plan_instruction, composite_mi, active_link_names, index);

          ompl::base::ScopedState<> start_state(sub_prob->simple_setup->getStateSpace());
          if (isJointWaypoint(start_waypoint) || isStateWaypoint(start_waypoint))
          {
            assert(checkJointPositionFormat(joint_names, start_waypoint));
            const Eigen::VectorXd& prev_position = getJointPosition(start_waypoint);
            cur_plan_profile->applyStartStates(
                *sub_prob, prev_position, *start_instruction, composite_mi, active_link_names, index);
          }
          else if (isCartesianWaypoint(start_waypoint))
          {
            const auto& prev_wp = start_waypoint.as<tesseract_planning::CartesianWaypoint>();
            cur_plan_profile->applyStartStates(
                *sub_prob, prev_wp, *start_instruction, composite_mi, active_link_names, index);
          }
          else
          {
            throw std::runtime_error("OMPLMotionPlannerDefaultConfig: unknown waypoint type");
          }

          problem.push_back(std::move(sub_prob));
          ++index;
        }
        else if (isCartesianWaypoint(plan_instruction.getWaypoint()))
        {
          const auto& cur_wp = plan_instruction.getWaypoint().as<tesseract_planning::CartesianWaypoint>();
          cur_plan_profile->applyGoalStates(
              *sub_prob, cur_wp, plan_instruction, composite_mi, active_link_names, index);

          if (index == 0)
          {
            ompl::base::ScopedState<> start_state(sub_prob->simple_setup->getStateSpace());
            if (isJointWaypoint(start_waypoint) || isStateWaypoint(start_waypoint))
            {
              assert(checkJointPositionFormat(joint_names, start_waypoint));
              const Eigen::VectorXd& prev_position = getJointPosition(start_waypoint);
              cur_plan_profile->applyStartStates(
                  *sub_prob, prev_position, *start_instruction, composite_mi, active_link_names, index);
            }
            else if (isCartesianWaypoint(start_waypoint))
            {
              const auto& prev_wp = start_waypoint.as<tesseract_planning::CartesianWaypoint>();
              cur_plan_profile->applyStartStates(
                  *sub_prob, prev_wp, *start_instruction, composite_mi, active_link_names, index);
            }
            else
            {
              throw std::runtime_error("OMPLMotionPlannerDefaultConfig: unknown waypoint type");
            }
          }
          else
          {
            /** @todo Update. Extract the solution for the previous plan and set as the start */
            assert(false);
          }

          problem.push_back(std::move(sub_prob));
          ++index;
        }
      }
      else
      {
        throw std::runtime_error("OMPLMotionPlannerDefaultConfig: Unsupported!");
      }

      start_waypoint = plan_instruction.getWaypoint(); /** @todo need to extract the solution for Cartesian waypoints
                                                           to work correctly*/
      start_instruction = &instruction;
    }
  }

  return problem;
}

}  // namespace tesseract_planning
