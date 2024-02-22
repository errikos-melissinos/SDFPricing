module ConsumptionFinance

#* import packages
import StochasticDiffEq as sde
import Statistics as stats
import Interpolations as intp
import DifferentialEquations.EnsembleAnalysis as ensemble
import Base.Iterators as iter

#* load zeroCouponSecurity
include("solveModel.jl")

export solveModel

end
