abstract type Solution end

struct SinglePayoffSolution{arr,intp,sp} <: Solution
    array::arr
    intp::intp
    samplePaths::sp
    problem::Problem
    solutionSettings::SolutionSettings
    outVariable::Int64

    SinglePayoffSolution(array, intp, samplePaths, problem, solutionSettings, outVariable) = new{typeof(array),typeof(intp),typeof(samplePaths)}(array, intp, samplePaths, problem, solutionSettings, outVariable)
end

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
