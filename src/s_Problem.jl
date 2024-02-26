struct Problem
    drift::Function
    diffusion::Function
    numNoiseVariables::Int64
    outVariables::Vector{Int64}
    terminalFunction::Function
    diagonalNoise::Bool

    #* Default Problem for testing purposes and convenience
    function Problem(;
        drift=((du, u, p, t) -> (du[1] = -0.1 * u[1]; du[2] = u[1])),
        diffusion=((du, u, p, t) -> (du[1] = 0.01; du[2] = 0.0)),
        numNoiseVariables=1,
        outVariables=[2],
        terminalFunction=(k, x, y, z) -> exp(-x),
        diagonalNoise=true
    )
        println("This is a default/example problem.")
        new(drift, diffusion, numNoiseVariables, outVariables, terminalFunction, diagonalNoise)
    end

    #* General constructor
    Problem(drift, diffusion, numNoiseVariables, outVariables, terminalFunction, diagonalNoise) = new(drift, diffusion, numNoiseVariables, outVariables, terminalFunction, diagonalNoise)

    #* adjust existing problem
    Problem(prob::Problem;
        drift=prob.drift,
        diffusion=prob.diffusion,
        numNoiseVariables=prob.numNoiseVariables,
        outVariables=prob.outVariables,
        terminalFunction=prob.terminalFunction,
        diagonalNoise=prob.diagonalNoise
    ) = new(drift, diffusion, numNoiseVariables, outVariables, terminalFunction, diagonalNoise)
end
