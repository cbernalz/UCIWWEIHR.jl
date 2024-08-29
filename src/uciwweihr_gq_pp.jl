# Fitting UCIWWEIHR model
# -------------------------------------------------
"""
    uciwweihr_gq_pp(...)
This generates quantities and a posterior predictive distribution for the bayesian semi-parametric model for the wastewater EIHR compartmental model.
The defaults for this fuction will follow those of the default simulation in generate_simulation_data_ww_eihr.jl function.

# Arguments
- `samples`: Samples from the posterior/prior distribution.
- `data_hosp`: An array of hospital data.
- `data_wastewater`: An array of pathogen genome concentration in localized wastewater data.  If this is not avaliable, the model used will be one that only uses hospital data.
- `obstimes`: An array of timepoints for observed hosp/wastewater.
- `param_change_times`: An array of timepoints where the parameters change.
- `seed::Int64=2024`: Seed for the random number generator.
- `params::uciwweihr_model_params`: A struct containing parameters for the model.
- `forecast::Bool=false`: A boolean to indicate if forecasting is to be done.
- `forecast_weeks::Int64=4`: Number of weeks to forecast.

# Returns
- Samples from the posterior or prior distribution.
"""

function uciwweihr_gq_pp(
    samples,
    data_hosp,
    data_wastewater;
    obstimes,
    param_change_times,
    seed::Int64=2024,
    params::uciwweihr_model_params,
    forecast::Bool=false, forecast_weeks::Int64=4
    )

    println("Using uciwweihr_model with wastewater!!!")
    obstimes = convert(Vector{Float64}, obstimes)
    param_change_times = convert(Vector{Float64}, param_change_times)


    if forecast
        last_value = obstimes[end]
        for i in 1:forecast_weeks
            next_value = last_value + 7
            push!(param_change_times, next_value)
            push!(obstimes, next_value)
            last_value = next_value
        end
        missing_data_ww = repeat([missing], length(obstimes))
        missing_data_hosp = repeat([missing], length(obstimes))
        data_hosp = vcat(data_hosp, repeat([data_hosp[end]], forecast_weeks))
        data_wastewater = vcat(data_wastewater, repeat([data_wastewater[end]], forecast_weeks))
    else
        missing_data_ww = repeat([missing], length(data_wastewater))
        missing_data_hosp = repeat([missing], length(data_hosp))
    end

    my_model = uciwweihr_model(
        data_hosp, 
        data_wastewater;
        obstimes, 
        param_change_times,
        params
    )


    #indices_to_keep = .!isnothing.(generated_quantities(my_model, samples))
    #samples_randn = ChainsCustomIndex(samples, indices_to_keep)

    #Random.seed!(seed)
    #gq_randn = Chains(generated_quantities(my_model, samples_randn))

    my_model_forecast_missing = uciwweihr_model(
        missing_data_hosp, 
        missing_data_ww;
        obstimes, 
        param_change_times,
        params
    )


    indices_to_keep = .!isnothing.(generated_quantities(my_model, samples))
    samples_randn = ChainsCustomIndex(samples, indices_to_keep)


    Random.seed!(seed)
    predictive_randn = predict(my_model_forecast_missing, samples_randn)

    Random.seed!(seed)
    gq_randn = Chains(generated_quantities(my_model, samples_randn))

    samples_df = DataFrame(samples)

    results = [DataFrame(predictive_randn), DataFrame(gq_randn), samples_df]


    return(results)
end

function uciwweihr_gq_pp(
    samples,
    data_hosp;
    obstimes,
    param_change_times,
    seed::Int64=2024,
    params::uciwweihr_model_params,
    forecast::Bool=false, forecast_weeks::Int64=4
    )
    println("Using uciwweihr_model without wastewater!!!")
    obstimes = convert(Vector{Float64}, obstimes)
    param_change_times = convert(Vector{Float64}, param_change_times)


    if forecast
        last_value = obstimes[end]
        for i in 1:forecast_weeks
            next_value = last_value + 7
            push!(param_change_times, next_value)
            push!(obstimes, next_value)
            last_value = next_value
        end
        missing_data_hosp = repeat([missing], length(obstimes))
        data_hosp = vcat(data_hosp, repeat([data_hosp[end]], forecast_weeks))
    else
        missing_data_hosp = repeat([missing], length(data_hosp))
    end

    my_model = uciwweihr_model(
        data_hosp;
        obstimes, 
        param_change_times,
        params
    )


    #indices_to_keep = .!isnothing.(generated_quantities(my_model, samples))
    #samples_randn = ChainsCustomIndex(samples, indices_to_keep)

    #Random.seed!(seed)
    #gq_randn = Chains(generated_quantities(my_model, samples_randn))

    my_model_forecast_missing = uciwweihr_model(
        missing_data_hosp;
        obstimes, 
        param_change_times,
        params
    )


    indices_to_keep = .!isnothing.(generated_quantities(my_model, samples))
    samples_randn = ChainsCustomIndex(samples, indices_to_keep)


    Random.seed!(seed)
    predictive_randn = predict(my_model_forecast_missing, samples_randn)

    Random.seed!(seed)
    gq_randn = Chains(generated_quantities(my_model, samples_randn))

    samples_df = DataFrame(samples)

    results = [DataFrame(predictive_randn), DataFrame(gq_randn), samples_df]


    return(results)
end