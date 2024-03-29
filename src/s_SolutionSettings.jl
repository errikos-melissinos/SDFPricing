
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

