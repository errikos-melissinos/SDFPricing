using Test
import ConsumptionFinance as cf
import DifferentialEquations.EnsembleAnalysis as ensemble
import StochasticDiffEq as sde
import Statistics as stats



paths = 5
xSpan = -0.01:0.01:0.01
tSpan = 0.0:0.25:1.0
(bondPrice,) = cf.zeroCouponSecurity(drift=(du, u, p, t) -> (du[1] = 0.1; du[2] = u[1]), diffusion=(du, u, p, t) -> (du[1] = 0.01; du[2] = 0.0), xSpans=(xSpan,), initialValues=[[x, 0.0] for x in xSpan], numNoiseVariables=1, outVariables=[2], terminalFunction=(ik, x, y, z) -> exp(-x), algorithm=sde.LambaEM(), pathsPerInitialValue=paths, tSpan=tSpan)
@test applicable(bondPrice, 0.0, 0.0)
@test length(bondPrice.examples) == length(xSpan)
@test size(bondPrice.solArray) == (length(tSpan), length(xSpan))
@test bondPrice(0.0, 0.0) ≈ 1.0

#* the expectation hypothesis approximately holds when volatility is very low
#* at the steady state yields are all equal.
r(x) = 0.01 + x #* short term interest rate function
volatility = 0.0000001
(bondPrice,) = cf.zeroCouponSecurity(drift=(du, u, p, t) -> (du[1] = -0.9u[1]; du[2] = r(u[1])), diffusion=(du, u, p, t) -> (du[1] = volatility; du[2] = 0.0), xSpans=(xSpan,), initialValues=[[x, 0.0] for x in xSpan], numNoiseVariables=1, outVariables=[2], terminalFunction=(ik, x, y, z) -> exp(-x), algorithm=sde.LambaEM(), pathsPerInitialValue=paths, tSpan=0.0:1.0:10.0)
bondYield(t, x) = -log(bondPrice(t, x)) / t
@test abs(bondYield(1.0, 0.0)) < r(0.0) + 0.00001
@test abs(bondYield(5.0, 0.0)) < r(0.0) + 0.00001
@test abs(bondYield(10.0, 0.0)) < r(0.0) + 0.00001

xSpans = (-0.01:0.01:0.01, -0.02:0.02:0.04)
tSpan = 0.0:1.0:10.0
y0 = 0.01
r2(x, y) = x + y
u0 = [[x, y, 0.0] for x in xSpans[1], y in xSpans[2]]
(bondPrice2,) = cf.zeroCouponSecurity(drift=(du, u, p, t) -> (du[1] = -0.92u[1]; du[2] = -0.88(u[2] - y0); du[3] = r2(u[1], u[2])), diffusion=(du, u, p, t) -> (du[1] = volatility; du[2] = volatility; du[3] = 0), xSpans=xSpans, initialValues=u0, numNoiseVariables=2, diagonalNoise=false, outVariables=[3], terminalFunction=(ik, x, y, z) -> exp(-x), algorithm=sde.LambaEM(), pathsPerInitialValue=paths, tSpan=tSpan)
@test applicable(bondPrice2, 0.0, 0.0, 0.0)
@test size(bondPrice2.examples) == length.(xSpans)
@test size(bondPrice2.solArray) == (length(tSpan), length.(xSpans)...)
@test bondPrice2(0.0, 0.0, 0.0) ≈ 1.0
@test bondPrice2(0.0, -0.01, 0.02) ≈ 1.0
bondYield2(t, x, y) = -log(bondPrice2(t, x, y)) / t
@test abs(bondYield2(1.0, 0.0, y0) - r2(0.0, y0)) < 0.00001
@test abs(bondYield2(5.0, 0.0, y0) - r2(0.0, y0)) < 0.00001
@test abs(bondYield2(10.0, 0.0, y0) - r2(0.0, y0)) < 0.00001

#* test initial values
@test size(bondPrice2.examples) == length.(xSpans)
@test bondPrice2.examples[1, 1].u[1][1] ≈ u0[1, 1][1]
@test bondPrice2.examples[1, 1].u[1][2] ≈ u0[1, 1][2]
@test bondPrice2.examples[1, 2].u[1][1] ≈ u0[1, 2][1]
@test bondPrice2.examples[1, 2].u[1][2] ≈ u0[1, 2][2]
