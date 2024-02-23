module ConsumptionFinance

#* import packages
import StochasticDiffEq as sde
import Statistics as stats
import Interpolations as intp
import DataInterpolations as dataIntp
import DifferentialEquations.EnsembleAnalysis as ensemble
import Base.Iterators as iter
import QuadGK

#* load basic solveModel
include("solveModel.jl")

#* load solvePerpetuity
include("solvePerpetuity.jl")


export solveModel, Model, solvePerpetuity

end
