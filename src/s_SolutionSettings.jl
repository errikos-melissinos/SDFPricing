"""
//////////////////////////////////////////////////////////////
struct SDFPricing.SolutionSettings{alg} struct
Variable that holds the information required for the solution
that are not strictly part of the problem.
==============================================================
Fields
======

xRanges                  :: Vector
        --> vector of ranges for the state variables. Just 
        one range in a vector, when there is only one state 
        variable.
initialValues            :: Vector{Vector{Float64}}
        --> vector of initial values for the state variables.
        There is a vector for each set of initial values,
        which contains the values, just one value if there is
        just one state variable.
tRange                   :: AbstractRange{Float64}
        --> time range from 0 up to the highest maturity of
        the bond for the SinglePayoffSolution. For the 
        ContinuousPayoffSolution, the time range should be 
        long so that the dividend strips reach a value of 0.
        And the price-dividend ratio is converged.
dt                       :: Float64
        --> time step for the simulation, it makes a 
        difference only when the algorithm required a fixed 
        time step. Otherwise, the algorithm adaptively picks 
        a time step.
pathsPerInitialValue     :: Int64
        --> number of paths to simulate for each initial 
        value.
algorithm                :: alg
        --> algorithm to use for the simulation according 
        to the DifferentialEquations package.
continuousPayoffRanges   :: Vector
        --> xRanges that will only apply to the continuous 
        payoff solution, by default equal to XRanges.
continuousPayoffVars     :: Vector
        --> specifies the variable for which the 
        price-dividend ratio is going to be computed. For 
        this variable, the SinglePayoffSolution is computed 
        as well, and then integrated to get the price-
        dividend ratio.
continuousPayoffDuration :: Float64
        --> specifies the period for which the dividend 
        strips are integrated. By default it is equal to 
        the end of the tRange. It can be specified to be 
        smaller if the user wants to calculate the price-
        dividend ratio for a security that pays for a 
        limited period of time.

* All arguments are given as keyword arguments. If no 
arguments are given, a default SolutionSettings variable
is created. 
==============================================================
==============================================================
"""
struct SolutionSettings{alg}
    xRanges::Vector
    initialValues::Vector{Vector{Float64}}
    tRange::AbstractRange{Float64}
    dt::Float64
    pathsPerInitialValue::Int64
    algorithm::alg
    continuousPayoffRanges::Vector
    continuousPayoffVars::Vector
    continuousPayoffDuration::Float64

    #* Default SolutionSettings for testing purposes and convenience and Constructor with keywords
    function SolutionSettings(;
        xRanges=[-0.03:0.01:0.05,],
        initialValues=[[x, 0.0] for x in -0.03:0.01:0.05],
        tRange=0.0:0.5:5.0,
        dt=0.5,
        algorithm=sde.EM(),
        pathsPerInitialValue=default_numPathsPerInitialValue,
        continuousPayoffRanges=xRanges,
        continuousPayoffVars=[],
        continuousPayoffDuration=tRange[end])

        new{typeof(algorithm)}(xRanges, initialValues, tRange, dt, pathsPerInitialValue, algorithm, continuousPayoffRanges, continuousPayoffVars, continuousPayoffDuration)
    end

    #* Classic constructor
    SolutionSettings(xRanges, initialValues, tRange, dt, numPathsPerInitialValue, algorithm, continuousPayoffRanges, continuousPayoffVars, continuousPayoffDuration) = new{typeof(algorithm)}(xRanges, initialValues, tRange, dt, numPathsPerInitialValue, algorithm, continuousPayoffRanges, continuousPayoffVars, continuousPayoffDuration)

    #* adjust existing settings
    SolutionSettings(sett::SolutionSettings;
        xRanges=sett.xRanges,
        initialValues=sett.initialValues,
        tRange=sett.tRange,
        dt=sett.dt,
        pathsPerInitialValue=sett.pathsPerInitialValue,
        algorithm=sett.algorithm,
        continuousPayoffRanges=sett.continuousPayoffRanges,
        continuousPayoffVars=sett.continuousPayoffVars,
        continuousPayoffDuration=sett.continuousPayoffDuration
    ) = new{typeof(sett.algorithm)}(xRanges, initialValues, tRange, dt, pathsPerInitialValue, algorithm, continuousPayoffRanges, continuousPayoffVars, continuousPayoffDuration)
end


