"""
//////////////////////////////////////////////////////////////
SDFPricing.Problem struct
Variable that holds the information of the problem.
==============================================================

Fields
======

drift             :: Function
        --> drift function of the SDE, consistent with 
        DifferentialEquations package.
diffusion         :: Function
        --> diffusion function of the SDE, consistent with 
        DifferentialEquations package.
numNoiseVariables :: Int64
        --> number of noise variables in the SDE
outVariables      :: Vector{Int64}
        --> vector of indices of the output variables 
        (corresponds to the simulation of the discounting 
        integral, in the drift and diffusion functions 
        of the problem struct).
terminalFunction  :: Function
        --> terminal function, which includes disounting. 
        In most cases default should be used unless, there 
        is a more exotic payoff at maturity (for example 
        depending on the state variable).
diagonalNoise     :: Bool
        --> true if the noise is diagonal (implies one noise 
        variable), should be false if more than one noise 
        variables are used.

* All arguments are given as keyword arguments. If no 
arguments are given, a default problem is created. 
==============================================================
==============================================================
"""
struct Problem
        drift::Function
        diffusion::Function
        numNoiseVariables::Int64
        outVariables::Vector{Int64}
        terminalFunction::Function
        diagonalNoise::Bool

        function Problem(;
                drift=((du, u, p, t) -> (du[1] = -0.1 * u[1]; du[2] = u[1])),
                diffusion=((du, u, p, t) -> (du[1] = 0.01; du[2] = 0.0)),
                numNoiseVariables=1,
                outVariables=[2],
                terminalFunction=(k, x, y, z) -> exp(-x),
                diagonalNoise=true
        )
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
