using Test
import SDFPricing as sdf
import DifferentialEquations.EnsembleAnalysis as ensemble
import StochasticDiffEq as sde
import Statistics as stats


paths = 5
xSpan = -0.01:0.01:0.01
tSpan = 0.0:0.25:1.0
prob = sdf.Problem(drift=(du, u, p, t) -> (du[1] = 0.1; du[2] = u[1]), diffusion=(du, u, p, t) -> (du[1] = 0.01; du[2] = 0.0), numNoiseVariables=1, outVariables=[2], terminalFunction=(ik, x, y, z) -> exp(-x))
sett = sdf.SolutionSettings(xRanges=[xSpan,], initialValues=[[x, 0.0] for x in xSpan], algorithm=sde.LambaEM(), pathsPerInitialValue=paths, tRange=tSpan)
((bondPrice,),) = sdf.solve(prob, sett)
@test applicable(bondPrice, 0.0, 0.0)
@test length(bondPrice.samplePaths) == length(xSpan)
@test size(bondPrice.array) == (length(tSpan), length(xSpan))
@test bondPrice(0.0, 0.0) ≈ 1.0

#* the expectation hypothesis approximately holds when volatility is very low
#* at the steady state yields are all equal.
r(x) = 0.01 + x #* short term interest rate function
volatility = 0.0000001
prob2 = sdf.Problem(prob; drift=(du, u, p, t) -> (du[1] = -0.9u[1]; du[2] = r(u[1])), diffusion=(du, u, p, t) -> (du[1] = volatility; du[2] = 0.0))
sett2 = sdf.SolutionSettings(sett; tRange=0.0:0.5:10.0)
((bondPrice2,),) = sdf.solve(prob2, sett2)
bondYield(t, x) = -log(bondPrice2(t, x)) / t
@test abs(bondYield(1.0, 0.0)) < r(0.0) + 0.00001
@test abs(bondYield(5.0, 0.0)) < r(0.0) + 0.00001
@test abs(bondYield(10.0, 0.0)) < r(0.0) + 0.00001


#* test two state variables
volatility = 0.0000001
xRanges = [-0.05:0.01:0.05, -0.05:0.01:0.04]
tRange = 0.0:1.0:10.0
y0 = 0.01
r3(x, y) = x + y
u0 = vcat([[x, y, 0.0] for y in xRanges[2] for x in xRanges[1]])
u0 = sdf.toVector(xRanges, [3])
prob3 = sdf.Problem(drift=(du, u, p, t) -> (du[1] = -0.08 * u[1]; du[2] = -0.12(u[2] - y0); du[3] = r3(u[1], u[2])), diffusion=(du, u, p, t) -> (du[1, 1] = volatility; du[1, 2] = 0.0; du[2, 1] = 0.0; du[2, 2] = volatility; du[3, 1] = 0; du[3, 2]), numNoiseVariables=2, outVariables=[3], diagonalNoise=false)
sett3 = sdf.SolutionSettings(xRanges=xRanges, initialValues=u0, tRange=tRange, algorithm=sde.LambaEM())
((bondPrice3,),) = sdf.solve(prob3, sett3)
@test applicable(bondPrice3, 0.0, 0.0, 0.0)
@test collect(size(bondPrice3.samplePaths)) == length.(xRanges)
@test size(bondPrice3.array) == (length(tRange), length.(xRanges)...)
@test bondPrice3(0.0, 0.0, 0.0) ≈ 1.0
@test bondPrice3(0.0, -0.01, 0.02) ≈ 1.0
bondYield3(t, x, y) = -log(bondPrice3(t, x, y)) / t
@test abs(bondYield3(1.0, 0.0, y0) - r3(0.0, y0)) < 0.00001
@test abs(bondYield3(5.0, 0.0, y0) - r3(0.0, y0)) < 0.00001
@test abs(bondYield3(10.0, 0.0, y0) - r3(0.0, y0)) < 0.00001
#* test initial values
@test collect(size(bondPrice3.samplePaths)) == length.(xRanges)
@test bondPrice3.samplePaths[1, 1][:, 1] ≈ u0[1]
@test bondPrice3.samplePaths[2, 1][:, 1] ≈ u0[2]
@test bondPrice3.samplePaths[1, 2][:, 1] ≈ u0[1+length(xRanges[1])]


#* test continuous payoff solution
sett4 = sdf.SolutionSettings(sett3; continuousPayoffVars=[3])
((bondPrice4,), (priceConsumptionRatio,)) = sdf.solve(prob3, sett4)

using Plots
x = collect(xRanges[1])
plot(x, priceConsumptionRatio.(x, 0.0))
