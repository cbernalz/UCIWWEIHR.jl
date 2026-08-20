# GQ & PP UCIWWEIHR model
# -------------------------------------------------
"""
    generate_pq_pp(...)
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
- `return_bool::Bool=true`: A boolean to indicate if the model is to use the return statement.  **Only set to false if only forecast is desired**
- `gq_bool::Bool=true`: A boolean to indicate if the model is to generate quantities.  **Only set to false if only forecast is desired**

# Returns
- Samples from the posterior or prior distribution.
"""
function generate_pq_pp(
    samples,
    data_hosp,
    data_wastewater,
    obstimes_hosp,
    obstimes_wastewater,
    param_change_times,
    params::model_params_time_var_hosp_prev;
    seed::Int64=2024,
    forecast::Bool=false, forecast_days::Int64=14,
    weekly_bool::Bool=false
)
## prevalence model
    println("Generating quantities using uciwweihr_model with wastewater and time-varying hospitalization probability - Prevalence Model!!!")
    obstimes_hosp = convert(Vector{Int64}, obstimes_hosp)
    obstimes_wastewater = convert(Vector{Int64}, obstimes_wastewater)
    param_change_times = convert(Vector{Int64}, param_change_times)
    obstimes = unique(vcat(obstimes_hosp, obstimes_wastewater))
    obstimes = sort(obstimes)
    
    if forecast
        last_value = obstimes_hosp[end]
        if weekly_bool
            obstimes_hosp = vcat(obstimes_hosp,(last_value+7):7:(last_value+forecast_days))
            obstimes_wastewater = vcat(obstimes_wastewater,(last_value+7):7:(last_value+forecast_days))
        else
            obstimes_hosp = vcat(obstimes_hosp,(last_value+1):(last_value+forecast_days))
            obstimes_wastewater = vcat(obstimes_wastewater,(last_value+1):(last_value+forecast_days))
        end
        missing_data_hosp = repeat([missing], length(obstimes_hosp))
        missing_data_ww = repeat([missing], length(obstimes_wastewater))
        data_hosp = vcat(data_hosp, repeat([data_hosp[end]], forecast_days))
        data_wastewater = vcat(data_wastewater, repeat([data_wastewater[end]], forecast_days))
    else
        missing_data_ww = repeat([missing], length(data_wastewater))
        missing_data_hosp = repeat([missing], length(data_hosp))
    end

    my_model = uciwweihr_model(
        data_hosp,
        data_wastewater,
        obstimes_hosp,
        obstimes_wastewater,
        param_change_times,
        params;
    )
    my_model_forecast_missing = uciwweihr_model(
        missing_data_hosp,
        missing_data_ww,
        obstimes_hosp,
        obstimes_wastewater,
        param_change_times,
        params;
    )

    samples_df = DataFrame(samples)

    indices_to_keep = .!isnothing.(generated_quantities(my_model, samples))
    samples_randn = ChainsCustomIndex(samples, indices_to_keep)

    Random.seed!(seed)
    predictive_randn = predict(my_model_forecast_missing, samples_randn)
    Random.seed!(seed)
    println("Generating quantities...")
    gq_randn = Chains(generated_quantities(my_model, samples_randn))
    results = [DataFrame(predictive_randn), DataFrame(gq_randn), samples_df]

    return(results)
    
end



function generate_pq_pp(
    samples,
    data_hosp,
    data_wastewater,
    obstimes_hosp,
    obstimes_wastewater,
    param_change_times,
    params::model_params_time_var_hosp_inc;
    seed::Int64=2024,
    forecast::Bool=false, forecast_days::Int64=14,
    weekly_bool::Bool=false
)
    ## incidence model
    println("Generating quantities using uciwweihr_model with wastewater and time-varying hospitalization probability - Incidence Model!!!")
    obstimes_hosp = convert(Vector{Int64}, obstimes_hosp)
    obstimes_wastewater = convert(Vector{Int64}, obstimes_wastewater)
    param_change_times = convert(Vector{Int64}, param_change_times)
    obstimes = unique(vcat(obstimes_hosp, obstimes_wastewater))
    obstimes = sort(obstimes)
    
    if forecast
        last_value = obstimes_hosp[end]
        if weekly_bool
            obstimes_hosp = vcat(obstimes_hosp,(last_value+7):7:(last_value+forecast_days))
            obstimes_wastewater = vcat(obstimes_wastewater,(last_value+7):7:(last_value+forecast_days))
        else
            obstimes_hosp = vcat(obstimes_hosp,(last_value+1):(last_value+forecast_days))
            obstimes_wastewater = vcat(obstimes_wastewater,(last_value+1):(last_value+forecast_days))
        end
        missing_data_hosp = repeat([missing], length(obstimes_hosp))
        missing_data_ww = repeat([missing], length(obstimes_wastewater))
        data_hosp = vcat(data_hosp, repeat([data_hosp[end]], forecast_days))
        data_wastewater = vcat(data_wastewater, repeat([data_wastewater[end]], forecast_days))
    else
        missing_data_ww = repeat([missing], length(data_wastewater))
        missing_data_hosp = repeat([missing], length(data_hosp))
    end

    my_model = uciwweihr_model(
        data_hosp,
        data_wastewater,
        obstimes_hosp,
        obstimes_wastewater,
        param_change_times,
        params;
    )
    my_model_forecast_missing = uciwweihr_model(
        missing_data_hosp,
        missing_data_ww,
        obstimes_hosp,
        obstimes_wastewater,
        param_change_times,
        params;
    )

    samples_df = DataFrame(samples)

    indices_to_keep = .!isnothing.(generated_quantities(my_model, samples))
    samples_randn = ChainsCustomIndex(samples, indices_to_keep)

    Random.seed!(seed)
    predictive_randn = predict(my_model_forecast_missing, samples_randn)
    Random.seed!(seed)
    println("Generating quantities...")
    gq_randn = Chains(generated_quantities(my_model, samples_randn))
    results = [DataFrame(predictive_randn), DataFrame(gq_randn), samples_df]

    return(results)
    
end

function generate_pq_pp(
    samples,
    data_hosp_prev,
    data_hosp_inc,
    data_wastewater,
    obstimes_hosp_prev,
    obstimes_hosp_inc,
    obstimes_wastewater,
    param_change_times,
    params::model_params_hosp_inc_prev;
    seed::Int64=2024,
    forecast::Bool=false, forecast_days::Int64=14,
    weekly_bool::Bool=false
)
## prevalence model
    println("Generating quantities using uciwweihr_model with wastewater and time-varying hospitalization probability - Incidence and Prevalence Model!!!")
    obstimes_hosp_prev = convert(Vector{Int64}, obstimes_hosp_prev)
    obstimes_hosp_inc = convert(Vector{Int64}, obstimes_hosp_inc)
    obstimes_wastewater = convert(Vector{Int64}, obstimes_wastewater)
    param_change_times = convert(Vector{Int64}, param_change_times)
    obstimes = unique(vcat(obstimes_hosp_prev, obstimes_hosp_inc, obstimes_wastewater))
    obstimes = sort(obstimes)
    
    if forecast
        last_value_prev = obstimes_hosp_prev[end]
        last_value_inc = obstimes_hosp_inc[end]
        if weekly_bool
            obstimes_hosp_prev = vcat(obstimes_hosp_prev,(last_value+7):7:(last_value+forecast_days))
            obstimes_hosp_inc = vcat(obstimes_hosp_inc,(last_value+7):7:(last_value+forecast_days))
            obstimes_wastewater = vcat(obstimes_wastewater,(last_value+7):7:(last_value+forecast_days))
        else
            obstimes_hosp_prev = vcat(obstimes_hosp_prev,(last_value+1):(last_value+forecast_days))
            obstimes_hosp_inc = vcat(obstimes_hosp_inc,(last_value+1):(last_value+forecast_days))
            obstimes_wastewater = vcat(obstimes_wastewater,(last_value+1):(last_value+forecast_days))
        end
        missing_data_hosp_prev = repeat([missing], length(obstimes_hosp_prev))
        missing_data_hosp_inc = repeat([missing], length(obstimes_hosp_inc))
        missing_data_ww = repeat([missing], length(obstimes_wastewater))
        data_hosp_prev = vcat(data_hosp_prev, repeat([data_hosp_prev[end]], forecast_days))
        data_hosp_inc = vcat(data_hosp_inc, repeat([data_hosp_inc[end]], forecast_days))
        data_wastewater = vcat(data_wastewater, repeat([data_wastewater[end]], forecast_days))
    else
        missing_data_ww = repeat([missing], length(data_wastewater))
        missing_data_hosp_prev = repeat([missing], length(data_hosp_prev))
        missing_data_hosp_inc = repeat([missing], length(data_hosp_inc))
    end

    my_model = uciwweihr_model(
        data_hosp_prev,
        data_hosp_inc,
        data_wastewater,
        obstimes_hosp_prev,
        obstimes_hosp_inc,
        obstimes_wastewater,
        param_change_times,
        params;
    )
    my_model_forecast_missing = uciwweihr_model(
        missing_data_hosp_prev,
        missing_data_hosp_inc,
        missing_data_ww,
        obstimes_hosp_prev,
        obstimes_hosp_inc,
        obstimes_wastewater,
        param_change_times,
        params;
    )

    samples_df = DataFrame(samples)

    indices_to_keep = .!isnothing.(generated_quantities(my_model, samples))
    samples_randn = ChainsCustomIndex(samples, indices_to_keep)

    Random.seed!(seed)
    predictive_randn = predict(my_model_forecast_missing, samples_randn)
    Random.seed!(seed)
    println("Generating quantities...")
    gq_randn = Chains(generated_quantities(my_model, samples_randn))
    results = [DataFrame(predictive_randn), DataFrame(gq_randn), samples_df]

    return(results)
    
end


function generate_pq_pp(
    samples,
    data_hosp,
    obstimes_hosp,
    param_change_times,
    params::model_params_time_var_hosp_inc_no_ww;
    seed::Int64=2024,
    forecast::Bool=false, forecast_days::Int64=14,
    weekly_bool::Bool=false
)
    ## incidence model
    println("Generating quantities using uciwweihr_model w/out wastewater and nontime-varying hospitalization probability - Incidence Model!!!")
    obstimes_hosp = convert(Vector{Int64}, obstimes_hosp)
    param_change_times = convert(Vector{Int64}, param_change_times)
    obstimes = unique(vcat(obstimes_hosp))
    obstimes = sort(obstimes)
    
    if forecast
        last_value = obstimes_hosp[end]
        if weekly_bool
            obstimes_hosp = vcat(obstimes_hosp,(last_value+7):7:(last_value+forecast_days))
        else
            obstimes_hosp = vcat(obstimes_hosp,(last_value+1):(last_value+forecast_days))
        end
        missing_data_hosp = repeat([missing], length(obstimes_hosp))
        data_hosp = vcat(data_hosp, repeat([data_hosp[end]], forecast_days))
    else
        missing_data_hosp = repeat([missing], length(data_hosp))
    end

    my_model = uciwweihr_model(
        data_hosp,
        obstimes_hosp,
        param_change_times,
        params;
    )
    my_model_forecast_missing = uciwweihr_model(
        missing_data_hosp,
        obstimes_hosp,
        param_change_times,
        params;
    )

    samples_df = DataFrame(samples)

    indices_to_keep = .!isnothing.(generated_quantities(my_model, samples))
    samples_randn = ChainsCustomIndex(samples, indices_to_keep)

    Random.seed!(seed)
    predictive_randn = predict(my_model_forecast_missing, samples_randn)
    Random.seed!(seed)
    println("Generating quantities...")
    gq_randn = Chains(generated_quantities(my_model, samples_randn))
    results = [DataFrame(predictive_randn), DataFrame(gq_randn), samples_df]

    return(results)
    
end



function generate_pq_pp(
    samples,
    data_hosp,
    obstimes_hosp,
    param_change_times,
    params::model_params_non_time_var_hosp_no_ww;
    seed::Int64=2024,
    forecast::Bool=false, forecast_days::Int64=14
)
    println("Generating quantities using uciwweihr_model w/out wastewater and non-time varying hospitalization probability!!!")
    obstimes_hosp = convert(Vector{Int64}, obstimes_hosp)
    param_change_times = convert(Vector{Int64}, param_change_times)
    
    if forecast
        last_value = obstimes_hosp[end]
        obstimes_hosp = vcat(obstimes_hosp,(last_value+1):(last_value+forecast_days))
        missing_data_hosp = repeat([missing], length(obstimes_hosp))
        data_hosp = vcat(data_hosp, repeat([data_hosp[end]], forecast_days))
    else
        missing_data_hosp = repeat([missing], length(data_hosp))
    end

    my_model = uciwweihr_model(
        data_hosp,
        obstimes_hosp,
        param_change_times,
        params;
    )
    my_model_forecast_missing = uciwweihr_model(
        missing_data_hosp,
        obstimes_hosp,
        param_change_times,
        params;
    )

    samples_df = DataFrame(samples)

    indices_to_keep = .!isnothing.(generated_quantities(my_model, samples))
    samples_randn = ChainsCustomIndex(samples, indices_to_keep)

    Random.seed!(seed)
    predictive_randn = predict(my_model_forecast_missing, samples_randn)
    Random.seed!(seed)
    println("Generating quantities...")
    gq_randn = Chains(generated_quantities(my_model, samples_randn))
    results = [DataFrame(predictive_randn), DataFrame(gq_randn), samples_df]

    return(results)
    
end

## model without wastewater and with time-varying hospitalization probability - not implemented
## model with wastewater and without time-varying hospitalization probability - not implemented