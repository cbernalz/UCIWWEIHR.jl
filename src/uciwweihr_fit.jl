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
- `seed::Int64=2024`: Seed for the random number generator.
- `params::uciwweihr_model_params`: A struct containing parameters for the model.

# Returns
- Samples from the posterior or prior distribution.
"""
function uciwweihr_fit(
    data_hosp,
    data_wastewater;
    obstimes,
    param_change_times,
    priors_only::Bool=false,
    n_samples::Int64=500, n_chains::Int64=1, seed::Int64=2024,
    params::uciwweihr_model_params
    )
    println("Using uciwweihr_model with wastewater!!!")
    obstimes = convert(Vector{Float64}, obstimes)
    param_change_times = convert(Vector{Float64}, param_change_times)


    my_model = uciwweihr_model(
        data_hosp, 
        data_wastewater;
        obstimes, 
        param_change_times,
        params
    )


    # Sample Posterior
    if priors_only
        Random.seed!(seed)
        samples = sample(my_model, Prior(), MCMCThreads(), 400, n_chains)
    else
        Random.seed!(seed)
        samples = sample(my_model, NUTS(), MCMCThreads(), n_samples, n_chains)
    end
    return(samples)
end

function uciwweihr_fit(
    data_hosp;
    obstimes,
    param_change_times,
    priors_only::Bool=false,
    n_samples::Int64=500, n_chains::Int64=1, seed::Int64=2024,
    params::uciwweihr_model_params
    )
    println("Using uciwweihr_model without wastewater!!!")
    obstimes = convert(Vector{Float64}, obstimes)
    param_change_times = convert(Vector{Float64}, param_change_times)


    my_model = uciwweihr_model(
        data_hosp;
        obstimes, 
        param_change_times,
        params
    )


    # Sample Posterior
    if priors_only
        Random.seed!(seed)
        samples = sample(my_model, Prior(), MCMCThreads(), 400, n_chains)
    else
        Random.seed!(seed)
        samples = sample(my_model, NUTS(), MCMCThreads(), n_samples, n_chains)
    end
    return(samples)
end