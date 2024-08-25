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
- `E_init_sd::Float64=50.0`: Standard deviation for the initial number of exposed individuals.
- `E_init_mean::Int64=200`: Mean for the initial number of exposed individuals.
- `I_init_sd::Float64=20.0`: Standard deviation for the initial number of infected individuals.
- `I_init_mean::Int64=100`: Mean for the initial number of infected individuals.
- `H_init_sd::Float64=5.0`: Standard deviation for the initial number of hospitalized individuals.
- `H_init_mean::Int64=20`: Mean for the initial number of hospitalized individuals.
- `gamma_sd::Float64=0.02`: Standard deviation for the rate of incubation.
- `log_gamma_mean::Float64=log(1/4)`: Mean for the rate of incubation on log scale.
- `nu_sd::Float64=0.02`: Standard deviation for the rate of leaving the infected compartment.
- `log_nu_mean::Float64=log(1/7)`: Mean for the rate of leaving the infected compartment on the log scale.
- `epsilon_sd::Float64=0.02`: Standard deviation for the rate of hospitalization recovery.
- `log_epsilon_mean::Float64=log(1/5)`: Mean for the rate of hospitalization recovery on the log scale.
- `rho_gene_sd::Float64=0.02`: Standard deviation for the rho prior.
- `log_rho_gene_mean::Float64=log(0.011)`: Mean for the row prior on log scale.
- `tau_sd::Float64=0.02`: Standard deviation for the scale/variation of the log scale data.
- `log_tau_mean::Float64=log(0.1)`: Mean for the scale/variation of the log scale data on log scale itself.
- `df_shape::Float64=2.0`: Shape parameter for the gamma distribution.
- `df_scale::Float64=10.0`: Scale parameter for the gamma distribution.
- `sigma_hosp_sd::Float64=50.0`: Standard deviation for the negative binomial distribution for hospital data.
- `sigma_hosp_mean::Float64=500.0`: Mean for the negative binomial distribution for hospital data.
- `Rt_init_sd::Float64=0.3`: Standard deviation for the initial value of the time-varying reproduction number.
- `Rt_init_mean::Float64=0.2`: Mean for the initial value of the time-varying reproduction number.
- `sigma_Rt_sd::Float64=0.2`: Standard deviation for normal prior of log time-varying reproduction number standard deviation.
- `sigma_Rt_mean::Float64=-3.0`: Mean for normal prior of log time-varying reproduction number standard deviation.
- `w_init_sd::Float64=0.1`: Standard deviation for the initial value of the time-varying hospitalization rate.
- `w_init_mean::Float64=log(0.35)`: Mean for the initial value of the time-varying hospitalization rate.
- `sigma_w_sd::Float64=0.2`: Standard deviation for normal prior of log time-varying hospitalization rate standard deviation.
- `sigma_w_mean::Float64=-3.5`: Mean for normal prior of time-varying hospitalization rate standard deviation.
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
    E_init_sd::Float64=50.0, E_init_mean::Int64=200,
    I_init_sd::Float64=20.0, I_init_mean::Int64=100,
    H_init_sd::Float64=5.0, H_init_mean::Int64=20,
    gamma_sd::Float64=0.02, log_gamma_mean::Float64=log(1/4),
    nu_sd::Float64=0.02, log_nu_mean::Float64=log(1/7),
    epsilon_sd::Float64=0.02, log_epsilon_mean::Float64=log(1/5),
    rho_gene_sd::Float64=0.02, log_rho_gene_mean::Float64=log(0.011),
    tau_sd::Float64=0.02, log_tau_mean::Float64=log(0.1),
    df_shape::Float64=2.0, df_scale::Float64=10.0,
    sigma_hosp_sd::Float64=50.0, sigma_hosp_mean::Float64=500.0,
    Rt_init_sd::Float64=0.3, Rt_init_mean::Float64=0.2,
    sigma_Rt_sd::Float64=0.2, sigma_Rt_mean::Float64=-3.0,
    w_init_sd::Float64=0.1, w_init_mean::Float64=log(0.35),
    sigma_w_sd::Float64=0.2, sigma_w_mean::Float64=-3.5,
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
        E_init_sd, E_init_mean,
        I_init_sd, I_init_mean,
        H_init_sd, H_init_mean,
        gamma_sd, log_gamma_mean,
        nu_sd, log_nu_mean,
        epsilon_sd, log_epsilon_mean,
        rho_gene_sd, log_rho_gene_mean,
        tau_sd, log_tau_mean,
        df_shape, df_scale,
        sigma_hosp_sd, sigma_hosp_mean,
        Rt_init_sd, Rt_init_mean,
        sigma_Rt_sd, sigma_Rt_mean,
        w_init_sd, w_init_mean,
        sigma_w_sd, sigma_w_mean
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
        E_init_sd, E_init_mean,
        I_init_sd, I_init_mean,
        H_init_sd, H_init_mean,
        gamma_sd, log_gamma_mean,
        nu_sd, log_nu_mean,
        epsilon_sd, log_epsilon_mean,
        rho_gene_sd, log_rho_gene_mean,
        tau_sd, log_tau_mean,
        df_shape, df_scale,
        sigma_hosp_sd, sigma_hosp_mean,
        Rt_init_sd, Rt_init_mean,
        sigma_Rt_sd, sigma_Rt_mean,
        w_init_sd, w_init_mean,
        sigma_w_sd, sigma_w_mean
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
    E_init_sd::Float64=50.0, E_init_mean::Int64=200,
    I_init_sd::Float64=20.0, I_init_mean::Int64=100,
    H_init_sd::Float64=5.0, H_init_mean::Int64=20,
    gamma_sd::Float64=0.02, log_gamma_mean::Float64=log(1/4),
    nu_sd::Float64=0.02, log_nu_mean::Float64=log(1/7),
    epsilon_sd::Float64=0.02, log_epsilon_mean::Float64=log(1/5),
    sigma_hosp_sd::Float64=50.0, sigma_hosp_mean::Float64=500.0,
    Rt_init_sd::Float64=0.3, Rt_init_mean::Float64=0.2,
    sigma_Rt_sd::Float64=0.2, sigma_Rt_mean::Float64=-3.0,
    w_init_sd::Float64=0.1, w_init_mean::Float64=log(0.35),
    sigma_w_sd::Float64=0.2, sigma_w_mean::Float64=-3.5,
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
        E_init_sd, E_init_mean,
        I_init_sd, I_init_mean,
        H_init_sd, H_init_mean,
        gamma_sd, log_gamma_mean,
        nu_sd, log_nu_mean,
        epsilon_sd, log_epsilon_mean,
        sigma_hosp_sd, sigma_hosp_mean,
        Rt_init_sd, Rt_init_mean,
        sigma_Rt_sd, sigma_Rt_mean,
        w_init_sd, w_init_mean,
        sigma_w_sd, sigma_w_mean
    )


    #indices_to_keep = .!isnothing.(generated_quantities(my_model, samples))
    #samples_randn = ChainsCustomIndex(samples, indices_to_keep)

    #Random.seed!(seed)
    #gq_randn = Chains(generated_quantities(my_model, samples_randn))

    my_model_forecast_missing = uciwweihr_model(
        missing_data_hosp;
        obstimes, 
        param_change_times,
        E_init_sd, E_init_mean,
        I_init_sd, I_init_mean,
        H_init_sd, H_init_mean,
        gamma_sd, log_gamma_mean,
        nu_sd, log_nu_mean,
        epsilon_sd, log_epsilon_mean,
        sigma_hosp_sd, sigma_hosp_mean,
        Rt_init_sd, Rt_init_mean,
        sigma_Rt_sd, sigma_Rt_mean,
        w_init_sd, w_init_mean,
        sigma_w_sd, sigma_w_mean
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