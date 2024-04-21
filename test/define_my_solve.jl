using Test
import SDFPricing as sdf
import DifferentialEquations.EnsembleAnalysis as ensemble
import StochasticDiffEq as sde
import Statistics as stats


#* test the defineMySolve function
paths, t0, t1, dt = 10, 0.0, 1.0, 0.05
alg = sde.EM()
tempFunc = sdf.defineMySolve(sde.EM(), dt, t0:2dt:t1)
W = sde.WienerProcess(0.0, 0.0, 0.0)
problem = sde.SDEProblem((du, u, p, t) -> (du[1] = 0.1), (du, u, p, t) -> (du[1] = 0.05), [0], (0.0, 1.0), noise=W)
ensembleProblem = ensemble.EnsembleProblem(problem)


sol = tempFunc(ensembleProblem, paths)
@test sol isa sde.EnsembleSolution
@test length(sol) == paths
@test sol[:, 1] isa sde.RODESolution
@test length(sol[:, 1]) ≈ t1 / 2dt + 1
@test applicable(sol[:, 1], 0.0)
@test typeof(sol[:, 1].alg) == typeof(alg)
@test sol[:, 1].t[1] == t0
@test sol[:, 1].t[end] == t1
@test sol[:, 1].t[2] - sol[:, 1].t[1] ≈ 2dt
@test stats.var(sol[:, 1].u)[1] > 0.00001

numInitialVals = 5
u0 = [[0.1 * i] for i in 1:numInitialVals]
problem = sde.SDEProblem((du, u, p, t) -> (du[1] = 0.1), (du, u, p, t) -> (du[1] = 0.01), [0], (0.0, 1.0), noise=W)
prob_func(prob, i, repeat) = sde.remake(prob, u0=u0[sdf.pickIndex(i, paths)])
ensembleProblem = ensemble.EnsembleProblem(problem, prob_func=prob_func)
sol = tempFunc(ensembleProblem, paths * numInitialVals)

#* test initial values
@test length(sol) == paths * numInitialVals
@test sol[:, 1].u[1][1] ≈ 0.1
@test sol[:, paths+1].u[1][1] ≈ 0.2
@test sol[:, 2paths+1].u[1][1] ≈ 0.3
@test sol[:, 3paths+1].u[1][1] ≈ 0.4
@test sol[:, 4paths+1].u[1][1] ≈ 0.5
@test sol[:, paths].u[1][1] ≈ 0.1
@test sol[:, 2paths].u[1][1] ≈ 0.2
@test sol[:, 3paths].u[1][1] ≈ 0.3
@test sol[:, 4paths].u[1][1] ≈ 0.4