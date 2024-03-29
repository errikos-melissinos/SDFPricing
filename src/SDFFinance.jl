module SDFFinance

#* import packages
import StochasticDiffEq as sde
import Statistics as stats
import Interpolations as intp
import DataInterpolations as dataIntp
import DifferentialEquations.EnsembleAnalysis as ensemble
import Base.Iterators as iter
import QuadGK

include("v_parameters.jl")
include("s_Problem.jl")
include("s_SolutionSettings.jl")
include("s_Solution.jl")
include("f_solve.jl")


export solve, Problem, SinglePayoffSolution, ContinuousPayoffSolution, SolutionSettings, toVector

end
