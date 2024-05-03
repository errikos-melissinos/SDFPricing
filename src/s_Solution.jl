abstract type Solution end

"""
//////////////////////////////////////////////////////////////
SDFPricing.SinglePayoffSolution{arr, intp, sp} struct
The variable returned by the solve function.
==============================================================

  Fields
  ======

  array            :: arr
        --> array containing prices for different states 
        and different maturities.
  intp             :: intp
        --> interpolation of prices as a function remaining 
        maturity and state, also called automatically when 
        the struct is called.
  samplePaths      :: sp
        --> a sample simulation.
  problem          :: SDFPricing.Problem
        --> the problem struct that was solved.
  solutionSettings :: SDFPricing.SolutionSettings
        --> the solutionSettings struct used to solve the
        problem.
  outVariable      :: Int64
        --> the index of the output variable (corresponds
        to the simulation of the discounting integral, in 
        the drift and diffusion functions of the problem 
        struct).

  Supertype Hierarchy
  ===================

  SDFPricing.SinglePayoffSolution{arr, intp, sp} 
  <: SDFPricing.Solution <: Any
  ==============================================================
  ==============================================================
  """
struct SinglePayoffSolution{arr,intp,sp} <: Solution
    array::arr
    intp::intp
    samplePaths::sp
    problem::Problem
    solutionSettings::SolutionSettings
    outVariable::Int64

    SinglePayoffSolution(array, intp, samplePaths, problem, solutionSettings, outVariable) = new{typeof(array),typeof(intp),typeof(samplePaths)}(array, intp, samplePaths, problem, solutionSettings, outVariable)
end

"""
//////////////////////////////////////////////////////////////
SDFPricing.ContinuousPayoffSolution{arr, intp} struct
The variable returned by the solve function.
==============================================================

  Fields
  ======

  array                :: arr
        --> array containing prices for different states and 
        different maturities
  intp                 :: intp
        --> interpolation of prices as a function remaining 
        maturity and state,
  singlePayoffSolution :: SDFPricing.SinglePayoffSolution
        --> contains the singlePayoffSolution struct that 
        is also computed during the solution
  problem              :: SDFPricing.Problem
        --> the problem struct that was solved
  solutionSettings     :: SDFPricing.SolutionSettings
        --> the solutionSettings struct used to solve the 
        problem
  outVariable          :: Int64
        --> the index of the output variable (corresponds to 
        the simulation of the discounting integral, in the 
        drift and diffusion functions of the problem struct)

  Supertype Hierarchy
  ===================

  SDFPricing.ContinuousPayoffSolution{arr, intp} 
  <: SDFPricing.Solution <: Any
  ==============================================================
  ==============================================================
"""
struct ContinuousPayoffSolution{arr,intp} <: Solution
    array::arr
    intp::intp
    singlePayoffSolution::SinglePayoffSolution
    problem::Problem
    solutionSettings::SolutionSettings
    outVariable::Int64

    ContinuousPayoffSolution(array, intp, singlePayoffSolution, problem, solutionSettings, outVariable) = new{typeof(array),typeof(intp)}(array, intp, singlePayoffSolution, problem, solutionSettings, outVariable)
end


#* add a call operator for the Problem struct that gives the price of the zero coupon security as a function of time and the state variable.
(sol::Solution)(args...) = sol.intp(args...)
