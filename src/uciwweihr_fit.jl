# Fitting UCIWWEIHR model
# -------------------------------------------------
"""
    uciwweihr_fit(...)
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
- `params::uciwweihr_model_params`: A struct containing parameters for the model.
- `init_params`: Initial parameters for the model.  If supplied, `uciwweihr_init_param` can supply these values.
- `return_bool`: A boolean to indicate if the model is to use the return statement.  **Only set to false if only forecast is desired**

# Returns
- Samples from the posterior or prior distribution.
"""
function uciwweihr_fit(
    data_hosp,
    data_wastewater,
    obstimes_hosp,
    obstimes_wastewater;
    param_change_times,
    priors_only::Bool=false,
    n_samples::Int64=500, n_chains::Int64=1, 
    n_discard_initial::Int64=0, seed::Int64=2024,
    params::uciwweihr_model_params,
    init_params = nothing,
    return_bool::Bool=true,
    )
    println("Using uciwweihr_model with wastewater!!!")
    obstimes_hosp = convert(Vector{Int64}, obstimes_hosp)
    obstimes_wastewater = convert(Vector{Int64}, obstimes_wastewater)
    param_change_times = convert(Vector{Int64}, param_change_times)

    my_model = uciwweihr_model(
        data_hosp, 
        data_wastewater,
        obstimes_hosp,
        obstimes_wastewater;
        param_change_times,
        params,
        return_bool
    )

    # Sample Posterior
    if priors_only
        Random.seed!(seed)
        samples = sample(my_model, Prior(), MCMCThreads(), 400, n_chains)
    else
        Random.seed!(seed)

        # Optimize
        if init_params === nothing
            samples = sample(my_model, NUTS(), MCMCThreads(), n_samples, n_chains, discard_initial = n_discard_initial, init_params = init_params)
        else
            samples = sample(my_model, NUTS(), MCMCThreads(), n_samples, n_chains, discard_initial = n_discard_initial)
        end
    end
    return(samples)
end

function uciwweihr_fit(
    data_hosp,
    obstimes_hosp;
    param_change_times,
    priors_only::Bool=false,
    n_samples::Int64=500, n_chains::Int64=1, 
    n_discard_initial::Int64=0, seed::Int64=2024,
    params::uciwweihr_model_params,
    init_params = nothing,
    return_bool::Bool=true,
    )
    println("Using uciwweihr_model without wastewater!!!")
    obstimes_hosp = convert(Vector{Int64}, obstimes_hosp)
    param_change_times = convert(Vector{Int64}, param_change_times)

    my_model_optimize = uciwweihr_model(
        data_hosp,
        obstimes_hosp; 
        param_change_times,
        params,
        return_bool
    )


    my_model = uciwweihr_model(
        data_hosp,
        obstimes_hosp; 
        param_change_times,
        params,
        return_bool
    )


    # Sample Posterior
    if priors_only
        Random.seed!(seed)
        samples = sample(my_model, Prior(), MCMCThreads(), 400, n_chains)
    else
        Random.seed!(seed)
        
        # Optimize
        if init_params === nothing
            samples = sample(my_model, NUTS(), MCMCThreads(), n_samples, n_chains, discard_initial = n_discard_initial, init_params = init_params)
        else
            samples = sample(my_model, NUTS(), MCMCThreads(), n_samples, n_chains, discard_initial = n_discard_initial)
        end
    end
    return(samples)
end