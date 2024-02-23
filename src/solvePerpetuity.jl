"""
define a struct for a model with a perpetuity.
"""
struct PerpetuityModel{Model,Intp}
    perpetuityIntp::Intp
    model::Model
    perpetuitySpans::Tuple
    perpetuityVar::Int64

    PerpetuityModel(intp, model, spans, var) = new{Model,Any}(intp, model, spans, var)
end

#* add a call operator for the perpetuityModel struct that gives the price consumption ratio of the perpetuity as a function the state variable.
(perpetuityModel::PerpetuityModel)(args...) = perpetuityModel.perpetuityIntp(args...)


function solvePerpetuity(; perpetuityVars=[2], perpetuitySpans=nothing, kwargs...)
    models = solveModel(; kwargs...)
    perpetuityModels = []

    if isnothing(perpetuitySpans)
        perpetuitySpans = get(kwargs, :xSpans, nothing)
    end
    for model in models
        if model.outVariable ∈ perpetuityVars
            integratedValues = [QuadGK.quadgk(t -> model(t, args...), 0.0, get(kwargs, :tSpan, nothing)[end])[1] for args in iter.product(perpetuitySpans...)]
            push!(perpetuityModels, PerpetuityModel(intp.CubicSplineInterpolation(perpetuitySpans..., integratedValues), model, perpetuitySpans, model.outVariable))
            # integratedValues = [QuadGK.quadgk(t -> model(t, args...), 0.0, get(kwargs, :tSpan, nothing)[end])[1] for args in iter.product(perpetuitySpans...)]
            # @show collect(iter.product(perpetuitySpans[1]))
            # @show integratedValues
            # @show perpetuitySpans
            # push!(perpetuityModels, PerpetuityModel(dataIntp.BSplineApprox(integratedValues, collect(iter.product(perpetuitySpans[1])), 3, 4, :Uniform, :Average, extrapolate=true), model, perpetuitySpans, model.outVariable))
        end
    end
    (perpetuityModels)
end
solvePerpetuity(model, perpetuityVars=[2]) = solvePerpetuity(drift=model.drift, diffusion=model.diffusion, xSpans=model.xSpans, initialValues=model.initialValues, numNoiseVariables=model.numNoiseVariables, outVariables=model.outVariables, terminalFunction=model.terminalFunction, diagonalNoise=model.diagonalNoise, tSpan=model.tSpan, pathsPerInitialValue=model.pathsPerInitialValue, dt=model.dt, algorithm=model.algorithm, perpetuityVars=perpetuityVars)
