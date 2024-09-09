"""
    repeated_forecast(...)
This is the function to make repreated forecast for a given forecast time span, `n_forecast_weeks`, and for given time points, `forecast_points`.
Plots can be made for these forecasts.  The output is an array of `uciwweihr_gq_pp` results for each `forecast_points`.

# Arguments
- `samples`: The MCMC samples from the model fit.
- `data_hosp`: The hospitalization data.
- `data_wastewater`: The wastewater data.
- `obstimes_hosp`: The time points for the hospitalization data.
- `obstimes_wastewater`: The time points for the wastewater data.
- `n_samples`: The number of samples to draw from the posterior.
- `param_change_times`: The time points where the parameters change.
- `params::uciwweihr_model_params`: The model parameters.
- `n_forecast_weeks`: The number of weeks to forecast.
- `forecast_points`: The time points to forecast, thees points should be present in obstimes_hosp.

# Returns
- An array of `uciwweihr_gq_pp` resuts and timeseries used for building for each `forecast_points`.
"""
function repeated_forecast(
    samples,
    data_hosp,
    data_wastewater,
    obstimes_hosp,
    obstimes_wastewater;
    n_samples::Int64,
    params::uciwweihr_model_params,
    n_forecast_weeks::Int64,
    forecast_points::Vector{Int64}
)
    results = []
    for max_point in forecast_points
        index_hosp = findfirst(x -> x == max_point, obstimes_hosp)
        if index_hosp === nothing
            error("THE FORECAST POINT SHOUDL BE PRESENT IN OBSTIMES_HOSP!!!")
        end
        index_ww = findfirst(x -> x == max_point, obstimes_wastewater)
        if index_ww === nothing
            index_ww = findfirst(x -> x < point, obstimes_wastewater)
            if index_ww === nothing
                error("FINDING THE INDEX FOR WW FORECAST POINT FAILED!!!")
            end
        end
        max_week = Int(ceil(max_point / 7))
        
        temp_data_hosp = data_hosp[1:index_hosp]
        temp_data_wastewater = data_wastewater[1:index_ww]
        temp_obstimes_hosp = obstimes_hosp[1:index_hosp]
        temp_obstimes_wastewater = obstimes_wastewater[1:index_ww]
        temp_param_change_times = 1:1:max_week
        temp_build_object = [
            temp_data_hosp,
            temp_data_wastewater,
            temp_obstimes_hosp,
            temp_obstimes_wastewater,
            temp_param_change_times
        ]
        
        samples = uciwweihr_fit(
            temp_data_hosp,
            temp_data_wastewater,
            temp_obstimes_hosp,
            temp_obstimes_wastewater;
            param_change_times = temp_param_change_times,
            priors_only = false,
            n_samples = n_samples,
            params = params
        )
        model_output = uciwweihr_gq_pp(
            samples,
            temp_data_hosp,
            temp_data_wastewater,
            temp_obstimes_hosp,
            temp_obstimes_wastewater;
            param_change_times = temp_param_change_times,
            params = params,
            forecast = true,
            forecast_weeks = n_forecast_weeks
        )
        push!(results, [temp_build_object, model_output])
    end
    return(results)
end
