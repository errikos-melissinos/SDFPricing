"""
This struct holds the info of the problem and (possibly) the results. It can also be used to get the price of the zero coupon security as a function of time and the state variable.
Its fields are:
- solArray:   the price of the security for the grid values of the state variable and time.
- solIntp:   the interpolation of the price of the security.
- drif
- diffusion
- xSpan
- initialValue
- numNoiseVariables::Int64
- outVariables::Vector{Int64}
- terminalFunctio
- diagonalNoise::Bool
- tSpa
- pathsPerInitialValue::Int64
- dt::Float64
- algorithm

This struct can also be initialised to not hold any results and then it can be passed to 
"""
struct Model{arr,intp,sol,dr,dif,xs,iv,tf,ts,alg}
    solArray::arr
    solIntp::intp
    examples::sol
    drift::dr
    diffusion::dif
    xSpans::xs
    initialValues::iv
    numNoiseVariables::Int64
    outVariables::Vector{Int64}
    terminalFunction::tf
    diagonalNoise::Bool
    tSpan::ts
    pathsPerInitialValue::Int64
    dt::Float64
    algorithm::alg

    #* Default constructor
    Model(;
        solArray=nothing,
        solIntp=nothing,
        examples=nothing,
        drift=((du, u, p, t) -> (du[1] = u[1]; du[2] = u[1])),
        diffusion=((du, u, p, t) -> (du[1] = 0.01; du[2] = 0.0)),
        xSpans=(-0.03:0.01:0.05,),
        initialValues=[[x, 0.0] for x in -0.03:0.01:0.05],
        tSpan=0.0:0.5:5.0,
        numNoiseVariables=1,
        outVariables=[2],
        terminalFunction=(k, x, y, z) -> exp(-x),
        dt=0.5,
        algorithm=sde.EM()
    ) = new{Any,Any,Any,Function,Function,Tuple,Array,Function,StepRangeLen,typeof(sde.EM())}(solArray, solIntp, examples, drift, diffusion, xSpans, initialValues, numNoiseVariables, outVariables, terminalFunction, true, tSpan, 8000, dt, algorithm)

    #* General constructor
    Model(means, interpolation, exampleSols, drift, diffusion, xSpans, initialValues, numNoiseVariables, outVariables, terminalFunction, diagonalNoise, tSpan, pathsPerInitialValue, dt, algorithm) = new{Any,Any,Any,Function,Function,Tuple,Array,Function,StepRangeLen,typeof(sde.EM())}(means, interpolation, exampleSols, drift, diffusion, xSpans, initialValues, numNoiseVariables, outVariables, terminalFunction, diagonalNoise, tSpan, pathsPerInitialValue, dt, algorithm)
end


#* add a call operator for the Model struct
function (model::Model)(t, args...)
    model.solIntp(t, args...)
end


#* Each initial value will be simulated pathsPerInitialValue times. This function will return the index of the initial value that will be used for the i-th simulation.
pickIndex(i, pathsPerInitialValue::Int64)::Int64 = 1 + div(i - 1, pathsPerInitialValue)

#* define the solve function generator to be used in the loops
#* the point is to predefine common arguments
defineMySolve(algorithm, dt, tSpan) = ((ensembleProblem, trajectories) -> sde.solve(ensembleProblem, algorithm, dt=dt, ensemble.EnsembleThreads(), saveat=tSpan, trajectories=trajectories))




"""
solveModel(drift,diffusion,xSpans,initialValues,numNoiseVariables::Int64,outVariables::Vector{Int64},terminalFunction,diagonalNoise::Bool=true,tSpan=0.0:1.0:20.0,pathsPerInitialValue::Int64=8000,dt::Float64=1 / 12, algorithm=sde.EM())

The pricing equation is:
E[d(ΛQ)/Λ]==0, where Λ is the SDF and Q is the price of the zero coupon security price. After applying Ito's lemma to dQ and given the SDF dΛ/Λ, we can derive a linear pde of the following form:
Q_m=-r(x)Q+drift(x)Q_x+diffusion(x)^2/2Q_xx
where Q_m is the derivative of Q with respect to matrutiy, Q_x is the derivative of Q with respect to x and Q_xx is the second derivative of Q with respect to x.
The solution to this pde is given by the Feynman-Kac formula
Q(m,x) = E[exp(-∫_0^m ỹ(x̃(s)) ds)*terminalFunction(x)|x(0)=x]
where dx̃ = drift(x̃)dt+diffusion(x̃)dW, x̃(0)=x 
and dỹ=r(x̃)dt, x̃(0)=x, ỹ(0)=0 and Q(t,x)=exp(-ỹ(t))
the terminalFunction is the payoff at maturity (typically equal to 1).
Applying Monte Carlo Simulations we can get the expectation shown above.

All arguments are entered with keywords:

drift: the drift function of the state variable(s), should follow the format of SDEProblem.
diffusion: the diffusion function of the state variable(s), should follow the format of SDEProblem.
xSpans: the range of the state variable(s), should be a tuple of ranges.
initialValues: the initial values of the state variable(s), should be a vector of vectors. 0.0 should also be included as an initial value for the auxiliary variable.
numNoiseVariables: the number of noise variables, should be 1 if the noise is diagonal.
outVariables: the indices of the variables that will be used in the terminal function.
terminalFunction: the terminal function needs to take four arguments, the first should be of type Val{k} where k is the index of the output variable, the second should be the corresponding state variable. If applicable the third is time and the fourth is the array of the all state variables. If not applicable just use some arguments as placeholders.
diagonalNoise: whether the noise is diagonal or not (true or false).
tSpan: the time span that is saved for the interpolation, this should be a range.
pathsPerInitialValue: the number of paths that will be simulated for each initial value.
dt: the time step that will be used in the simulation (this applies only for some algorithms).
algorithm: the algorithm that will be used in the simulation, the algorithms are from the StochasticDiffEq package.

Output
The output is a tuple that contains model structs. So if there is only one outVariable then the function should be called as follows:

(price,) = solveModel(...)
Price then can also be used as a function of remaining maturity and the state variable(s).
"""
function solveModel(;
    drift,
    diffusion,
    xSpans,
    initialValues,
    numNoiseVariables::Int64,
    outVariables::Vector{Int64},
    terminalFunction,
    diagonalNoise::Bool=true,
    tSpan=0.0:1.0:20.0,
    pathsPerInitialValue::Int64=8000,
    dt::Float64=1 / 12,
    algorithm=sde.EM())


    #* predefine the solve function that will be used several times
    # mySolve(ensembleProblem, trajectories) = sde.solve(ensembleProblem, algorithm, dt=dt, ensemble.EnsembleThreads(), saveat=tSpan, trajectories=trajectories)
    mySolve = defineMySolve(algorithm, dt, tSpan)

    #* define problem depending on whether the noise is diagonal or not
    if diagonalNoise
        W = sde.WienerProcess(0.0, 0.0, 0.0)
        problem = sde.SDEProblem(drift, diffusion, initialValues[1], (tSpan[1], tSpan[end]), noise=W)
    else
        if numNoiseVariables == 1
            println("You are using a non-diagonal noise process with only one noise variable. It is probably preferable to adjust this choice.")
        else
            problem = sde.SDEProblem(drift, diffusion, initialValues[1], (tSpan[1], tSpan[end]), noise_rate_prototype=zeros(length(initialValues[1]), numNoiseVariables))
        end
    end

    function prob_func(prob, i, repeat, initialValues, pathsPerInitialValue)
        prob = sde.remake(prob, u0=initialValues[pickIndex(i, pathsPerInitialValue)])
        prob
    end
    prob_func(prob, i, repeat) = prob_func(prob, i, repeat, initialValues, pathsPerInitialValue)

    # solve the problem
    ensembleProblem = ensemble.EnsembleProblem(problem, prob_func=prob_func)
    solution = mySolve(ensembleProblem, pathsPerInitialValue * length(initialValues))

    # get the FeynmanKacVariables
    values = Array{Float64}(undef, (length(tSpan), length(initialValues), pathsPerInitialValue))
    means = Array{Float64}(undef, (length(tSpan), length(initialValues)))
    models = []
    #* define myShape to make sure that the shape is always consistent
    myShape = (length(tSpan), length.(xSpans)...)
    for (ik, k) in enumerate(outVariables)
        for i in 1:pathsPerInitialValue
            for j in 1:length(initialValues)
                for t in 1:length(tSpan)
                    iPos = (j - 1) * pathsPerInitialValue + i
                    values[t, j, i] = terminalFunction(Val(ik), solution[iPos][k, t], solution[iPos].t[t], solution[iPos][:, t])
                end
            end
        end
        means .= stats.mean(values, dims=3)
        means2 = reshape(means, myShape...)
        exampleSols = reshape([solution[pathsPerInitialValue*i] for (i, _) in enumerate(initialValues)], myShape[2:end]...)
        temp = intp.scale(intp.interpolate(means2, intp.BSpline(intp.Cubic())), (tSpan, xSpans...))
        push!(models, Model(means2, temp, exampleSols, drift, diffusion, xSpans, initialValues, numNoiseVariables, outVariables, terminalFunction, diagonalNoise, tSpan, pathsPerInitialValue, dt, solution[1].alg))
    end

    models
end

solveModel(model) = solveModel(drift=model.drift, diffusion=model.diffusion, xSpans=model.xSpans, initialValues=model.initialValues, numNoiseVariables=model.numNoiseVariables, outVariables=model.outVariables, terminalFunction=model.terminalFunction, diagonalNoise=model.diagonalNoise, tSpan=model.tSpan, pathsPerInitialValue=model.pathsPerInitialValue, dt=model.dt, algorithm=model.algorithm)

