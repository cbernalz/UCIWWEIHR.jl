# Fitting UCIWWEIHR model
# -------------------------------------------------
"""
    fit(...)
This is the sampler for the bayesian semi-parametric model for the wastewater EIHR compartmental model.  
The defaults for this fuction will follow those of the default simulation in generate_simulation_data_ww_eihr.jl function.

# Arguments
- `data_hosp`: An array of hospital data.
- `data_wastewater`: An array of pathogen genome concentration in localized wastewater data.  If this is not avaliable, the model used will be one that only uses hospital data.
- `obstimes`: An array of timepoints for observed hosp/wastewater.
- `priors_only::Bool=false`: A boolean to indicate if only priors are to be sampled.
- `n_samples::Int64=500`: Number of samples to be drawn.
- `n_chains::Int64=1`: Number of chains to be run.
- `n_discard_initial::Int64=0`: Number of samples to be discarded.
- `seed::Int64=2024`: Seed for the random number generator.
- `params::uciwweihr_model_params#`: A struct containing parameters for the model.  # is either 1 or 2 for prior or hardcoded sigma_ww and sigma_hosp.
- `init_params`: Initial parameters for the model.  If supplied, `uciwweihr_init_param` can supply these values.
- `return_bool`: A boolean to indicate if the model is to use the return statement.  **Only set to false if only forecast is desired**

# Returns
- Samples from the posterior or prior distribution.
"""
function fit(
    data_hosp,
    data_wastewater,
    obstimes_hosp,
    obstimes_wastewater,
    param_change_times,
    params::model_params_time_var_hosp_prev;
    priors_only::Bool=false,
    n_samples::Int64=500, n_chains::Int64=1,
    n_discard_initial::Int64=0, seed::Int64=2024,
    init_params=nothing
    )
    ## prevalence model
    println("Fitting using uciwweihr_model with wastewater - Prevalence Model!!!")
    obstimes_hosp = convert(Vector{Int64}, obstimes_hosp)
    obstimes_wastewater = convert(Vector{Int64}, obstimes_wastewater)
    param_change_times = convert(Vector{Int64}, param_change_times)
    param_change_times = vcat(0, param_change_times)
    my_model = uciwweihr_model(
        data_hosp,
        data_wastewater,
        obstimes_hosp,
        obstimes_wastewater,
        param_change_times,
        params;
    )
    # Sample Posterior
    if priors_only
        Random.seed!(seed)
        samples = sample(my_model, Prior(), MCMCThreads(), 400, n_chains)
    else
        Random.seed!(seed)
        # Optimize
        if init_params === nothing
            samples = sample(my_model, NUTS(), MCMCThreads(), n_samples, n_chains, discard_initial = n_discard_initial)
        else
            println("Using Initial Parameters...")
            samples = sample(my_model, NUTS(), MCMCThreads(), n_samples, n_chains, discard_initial = n_discard_initial, init_params = init_params)
        end
    end 
    return(samples)
 
end


function fit(
    data_hosp,
    data_wastewater,
    obstimes_hosp,
    obstimes_wastewater,
    param_change_times,
    params::model_params_time_var_hosp_inc;
    incidence_model_bool=false,
    priors_only::Bool=false,
    n_samples::Int64=500, n_chains::Int64=1,
    n_discard_initial::Int64=0, seed::Int64=2024,
    init_params=nothing
    )
    ## incidence model
    println("Fitting using uciwweihr_model with wastewater - Incidence Model!!!")
    obstimes_hosp = convert(Vector{Int64}, obstimes_hosp)
    obstimes_wastewater = convert(Vector{Int64}, obstimes_wastewater)
    param_change_times = convert(Vector{Int64}, param_change_times)
    param_change_times = vcat(0, param_change_times)
    my_model = uciwweihr_model(
        data_hosp,
        data_wastewater,
        obstimes_hosp,
        obstimes_wastewater,
        param_change_times,
        params;
    )
    # Sample Posterior
    if priors_only
        Random.seed!(seed)
        samples = sample(my_model, Prior(), MCMCThreads(), 400, n_chains)
    else
        Random.seed!(seed)
        # Optimize
        if init_params === nothing
            samples = sample(my_model, NUTS(), MCMCThreads(), n_samples, n_chains, discard_initial = n_discard_initial)
        else
            println("Using Initial Parameters...")
            samples = sample(my_model, NUTS(), MCMCThreads(), n_samples, n_chains, discard_initial = n_discard_initial, init_params = init_params)
        end
    end 
    return(samples)
 
end



function fit(
    data_hosp,
    obstimes_hosp,
    param_change_times,
    params::model_params_non_time_var_hosp_no_ww;
    priors_only::Bool=false,
    n_samples::Int64=500, n_chains::Int64=1,
    n_discard_initial::Int64=0, seed::Int64=2024,
    init_params=nothing
    )
    println("Fitting using uciwweihr_model w/out wastewater.  With time-varying hospitalization probability!!!")
    obstimes_hosp = convert(Vector{Int64}, obstimes_hosp)
    param_change_times = convert(Vector{Int64}, param_change_times)
    param_change_times = vcat(0, param_change_times)
    obstimes = unique(vcat(obstimes_hosp))
    obstimes = sort(obstimes)
    my_model = uciwweihr_model(
        data_hosp,
        obstimes_hosp,
        param_change_times,
        params;
    )
    # Sample Posterior
    if priors_only
        Random.seed!(seed)
        samples = sample(my_model, Prior(), MCMCThreads(), 400, n_chains)
    else
        Random.seed!(seed)
        # Optimize
        if init_params === nothing
            samples = sample(my_model, NUTS(), MCMCThreads(), n_samples, n_chains, discard_initial = n_discard_initial)
        else
            println("Using Initial Parameters...")
            samples = sample(my_model, NUTS(), MCMCThreads(), n_samples, n_chains, discard_initial = n_discard_initial, init_params = init_params)
        end
    end 
    return(samples)
 
end

## model w/out wastewater and with non time-varying hospitalization probability - not implemented
## model w wastewater and with non time-varying hospitalization probability - not implemented