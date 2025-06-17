"""
Each initial value will be simulated pathsPerInitialValue times. 
This function will return the index of the initial value that will 
    be used for the i-th simulation.
"""
pickIndex(i, pathsPerInitialValue::Int64)::Int64 = 1 + div(i - 1, pathsPerInitialValue)

"""
Define the solve function generator to be used in the loops
the point is to predefine common arguments
"""
defineMySolve(algorithm, dt, tRange) = ((ensembleProblem, trajectories) -> sde.solve(ensembleProblem, algorithm, dt=dt, ensemble.EnsembleThreads(), saveat=tRange, trajectories=trajectories))

"""
Convenience function to turn ranges into arrays of values with the correct order
"""
function toVector(ranges, zeroPositions=[])
    function addZeros(values, zeroPositions)
        for pos in zeroPositions
            insert!(values, pos, 0.0)
        end
        values
    end
    reshape([addZeros(collect(values), zeroPositions) for values in Base.Iterators.product(ranges...)], *(length.(ranges)...))
end

"""
//////////////////////////////////////////////////////////////
SDFPricing.solve function
Main function that solves the pricing problem.
==============================================================

Input
=====
1. Problem struct 
2. SolutionSettings struct.

Output
======
Tuple of size two containing a:
1. tuple for each SinglePayoffSolution struct
2. tuple for each ContinuousPayoffSolution struct
* examples:
- if there is only one bond computed:
((bondPrice,),) = sdf.solve(problem, solutionSettings) 
- if there is one bond and one price-dividend ratio:
((bondPrice,dividendStrip),(pdratio,)) = sdf.solve(problem, solutionSettings) 
#* the dividend strip is needs to be computed as 
#* SinglePayoffSolution for the price-dividend ratio to also 
#* get computed.

* For convenience it is also possible to call solve without
any arguments, in which case a default solution is computed.
==============================================================
==============================================================
"""
solve() = solve(Problem(), SolutionSettings())

function solve(problem, solutionSettings)
    if length(solutionSettings.continuousPayoffRanges) == 0
        return (solve_singlePayoff(problem, solutionSettings),)
    else
        return solve_continuousPayoff(solve_singlePayoff(problem, solutionSettings))
    end
end

"""
Sub-function to solve the continuous payoff problem
"""
function solve_continuousPayoff(solutions)
    continuousPayoffSolutions = []
    for sol in solutions
        if sol.outVariable ∈ sol.solutionSettings.continuousPayoffVars
            integratedValues = [QuadGK.quadgk(t -> sol(t, args...), 0.0, sol.solutionSettings.continuousPayoffDuration)[1] for args in iter.product(sol.solutionSettings.continuousPayoffRanges...)]
            if length(sol.solutionSettings.continuousPayoffRanges) == 1
                push!(continuousPayoffSolutions, ContinuousPayoffSolution(integratedValues, dataIntp.BSplineApprox(integratedValues,
                        sol.solutionSettings.continuousPayoffRanges[1], 3, 4, :Uniform, :Average, extrapolation=dataIntp.ExtrapolationType.Constant), sol, sol.problem, sol.solutionSettings, sol.outVariable))
            else
                push!(continuousPayoffSolutions, ContinuousPayoffSolution(integratedValues, intp.cubic_spline_interpolation(Tuple(sol.solutionSettings.continuousPayoffRanges), integratedValues), sol, sol.problem, sol.solutionSettings, sol.outVariable))
            end

        end
    end
    (solutions, continuousPayoffSolutions)
end

"""
Internal sub-function to solve the single payoff problem
"""
function solve_singlePayoff(problem, solutionSettings)

    #* get the values of the inputs
    (drift, diffusion, numNoiseVariables, outVariables, terminalFunction, diagonalNoise) = (problem.drift, problem.diffusion, problem.numNoiseVariables, problem.outVariables, problem.terminalFunction, problem.diagonalNoise)
    (xRanges, initialValues, tRange, pathsPerInitialValue, dt, algorithm) = (solutionSettings.xRanges, solutionSettings.initialValues, solutionSettings.tRange, solutionSettings.pathsPerInitialValue, solutionSettings.dt, solutionSettings.algorithm)

    @assert *(length.(xRanges)...) == length(initialValues) "When solving the single payoff problem the xRanges should correspond to the same initialValues."


    #* predefine the solve function that will be used several times
    # mySolve(ensembleProblem, trajectories) = sde.solve(ensembleProblem, algorithm, dt=dt, ensemble.EnsembleThreads(), saveat=tRange, trajectories=trajectories)
    mySolve = defineMySolve(algorithm, dt, tRange)

    #* define problem depending on whether the noise is diagonal or not
    if diagonalNoise
        W = sde.WienerProcess(0.0, 0.0, 0.0)
        sdeProblem = sde.SDEProblem(drift, diffusion, initialValues[1], (tRange[1], tRange[end]), noise=W)
    else
        if numNoiseVariables == 1
            println("You are using a non-diagonal noise process with only one noise variable. It is probably preferable to adjust this choice.")
        else
            sdeProblem = sde.SDEProblem(drift, diffusion, initialValues[1], (tRange[1], tRange[end]), noise_rate_prototype=zeros(length(initialValues[1]), numNoiseVariables))
        end
    end

    function prob_func(prob, i, repeat, initialValues, pathsPerInitialValue)
        prob = sde.remake(prob, u0=initialValues[pickIndex(i, pathsPerInitialValue)])
        prob
    end
    prob_func(prob, i, repeat) = prob_func(prob, i, repeat, initialValues, pathsPerInitialValue)

    # solve the problem
    ensembleProblem = ensemble.EnsembleProblem(sdeProblem, prob_func=prob_func)
    sol = mySolve(ensembleProblem, pathsPerInitialValue * length(initialValues))

    # get the FeynmanKacVariables
    myShape = (length(tRange), length.(xRanges)...)
    means = Array{Float64}(undef, myShape)

    solutions = []
    indices = toVector((length.(xRanges) .|> x -> 1:x))
    for k in outVariables
        j = 0
        for ind in indices
            j += 1
            for t in 1:length(tRange)
                means[t, ind...] = solve_takeExpectation(pathsPerInitialValue, sol, k, t, j, terminalFunction, ind...)
            end
        end
        samples = reshape([sol[:, pathsPerInitialValue*i] for (i, _) in enumerate(initialValues)], myShape[2:end]...)
        interpolation = intp.scale(intp.interpolate(means, intp.BSpline(intp.Cubic())), (tRange, xRanges...))
        push!(solutions, SinglePayoffSolution(means, interpolation, samples, problem, solutionSettings, k))
    end
    solutions
end

"""
Internal sub-function to take hhe average corresponding to the Feynman-Kac formula
"""
@inline function solve_takeExpectation(pathsPerInitialValue, solution, k, t, j, terminalFunction, indices...)
    value = 0.0
    for (_, iPos) in enumerate(((j-1)*pathsPerInitialValue+1):(j*pathsPerInitialValue))
        value += terminalFunction(Val(k), solution[:, iPos].u[t][k], solution[:, iPos].t[t], solution[:, iPos].u[t])
    end
    value / pathsPerInitialValue
end

"""
//////////////////////////////////////////////////////////////
SDFPricing.derivatives function
Computes derivatives of price-dividend interpolation.
==============================================================

The input it the ContinuousPayoffSolution struct.
The output is a tuple of two functions,
the first and second derivatives.

The first derivative is computed with the built-in function 
from the DataInterpolations package.
The second derivative is then computed ``manually'' by 
computing a ratio of small differeneces.

The epsilon value can be chosen by the user as a keyword 
argument.
By default: epsilon=1e-5.
==============================================================
==============================================================
"""
function derivatives(sol; epsilon=1e-5)
    DPC(x) = dataIntp.derivative(sol.intp, x)
    D2PC(x) = x > (sol.solutionSettings.xRanges[1][end] - sol.solutionSettings.xRanges[1][1]) / 2 ? (-DPC(x - epsilon) + DPC(x)) / epsilon : (DPC(x + epsilon) - DPC(x)) / epsilon

    return (DPC, D2PC)
end

"""
Convenience function to simulate the state variable.

The inputs are:
- the solution struct (SinglePayoffSolution or ContinuousPayoffSolution)
- the time range
- the time step
- the initial values of the state variable(s)
"""
function simulate(sol, tRange, dt, initialValues)
    #* get the values of the inputs
    (drift, diffusion, diagonalNoise, algorithm, numNoiseVariables) = (sol.problem.drift, sol.problem.diffusion, sol.problem.diagonalNoise, sol.solutionSettings.algorithm, sol.problem.numNoiseVariables)

    #* define problem depending on whether the noise is diagonal or not
    if diagonalNoise
        W = sde.WienerProcess(0.0, 0.0, 0.0)
        sdeProblem = sde.SDEProblem(drift, diffusion, initialValues[1], (tRange[1], tRange[end]), noise=W)
    else
        if numNoiseVariables == 1
            println("You are using a non-diagonal noise process with only one noise variable. It is probably preferable to adjust this choice.")
        else
            sdeProblem = sde.SDEProblem(drift, diffusion, initialValues[1], (tRange[1], tRange[end]), noise_rate_prototype=zeros(length(initialValues[1]), numNoiseVariables))
        end
    end

    sde.solve(sdeProblem, algorithm, dt=dt, saveat=dt)

end
