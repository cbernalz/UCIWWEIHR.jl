# Fitting UCIWWEIHR model
# -------------------------------------------------
"""
    uciwweihr_init_param(...)
Gets initial parameters values for the UCIWWEIHR model.  Only need to run once.

# Arguments
- `data_hosp`: An array of hospital data.
- `data_wastewater`: An array of pathogen genome concentration in localized wastewater data.  If this is not avaliable, the model used will be one that only uses hospital data.
- `obstimes_hosp`: An array of timepoints for observed hosp.
- `obstimes_wastewater`: An array of timepoints for observed wastewater.
- `param_change_times`: An array of timepoints where the parameters change.
- `n_chains::Int64=1`: Number of chains to be run.
- `seed::Int64=2024`: Seed for the random number generator.
- `params::uciwweihr_model_params`: A struct containing parameters for the model.
- `verbose_optimize::Bool=false`: A boolean to indicate if the optimization process is to be verbose.

# Returns
- Samples from the posterior or prior distribution.
"""
function uciwweihr_init_param(
    data_hosp,
    data_wastewater,
    obstimes_hosp,
    obstimes_wastewater,
    param_change_times,
    params::uciwweihr_model_params2;
    n_chains::Int64=1,
    seed::Int64=2024,
    verbose_optimize::Bool=false
)
    println("Getting init param values...")
    println("Using uciwweihr_model with wastewater.  Priors on sigma_ww and sigma_hosp!!!")
    obstimes_hosp = convert(Vector{Int64}, obstimes_hosp)
    obstimes_wastewater = convert(Vector{Int64}, obstimes_wastewater)
    param_change_times = convert(Vector{Int64}, param_change_times)

    my_model_optimize = uciwweihr_model(
        data_hosp,
        data_wastewater,
        obstimes_hosp,
        obstimes_wastewater,
        obstimes,
        param_change_times,
        params;
    )
    Random.seed!(seed)
    # Optimize
    MAP_init = optimize_many_MAP2(my_model_optimize, 10, 1, verbose_optimize)[1]
    Random.seed!(seed)
    MAP_noise = vcat(randn(length(MAP_init) - 1, n_chains), transpose(zeros(n_chains)))
    MAP_noise = [MAP_noise[:,i] for i in 1:size(MAP_noise,2)]
    init = repeat([MAP_init], n_chains) .+ 0.05 * MAP_noise
    return(init)
end



function uciwweihr_init_param(
    data_hosp,
    obstimes_hosp,
    param_change_times,
    params::uciwweihr_model_params4;
    n_chains::Int64=1,
    seed::Int64=2024,
    verbose_optimize::Bool=false
)
    println("Getting init param values...")
    println("Using uciwweihr_model w/ wastewater.  Priors on sigma_hosp and const hospitalization probability!!!")
    obstimes_hosp = convert(Vector{Int64}, obstimes_hosp)
    param_change_times = convert(Vector{Int64}, param_change_times)

    my_model_optimize = uciwweihr_model(
        data_hosp,
        obstimes_hosp,
        param_change_times,
        params;
    )
    Random.seed!(seed)
    # Optimize
    MAP_init = optimize_many_MAP2(my_model_optimize, 10, 1, verbose_optimize)[1]
    Random.seed!(seed)
    MAP_noise = vcat(randn(length(MAP_init) - 1, n_chains), transpose(zeros(n_chains)))
    MAP_noise = [MAP_noise[:,i] for i in 1:size(MAP_noise,2)]
    init = repeat([MAP_init], n_chains) .+ 0.05 * MAP_noise
    return(init)
end


function uciwweihr_init_param(
    data_hosp,
    data_wastewater,
    obstimes_hosp,
    obstimes_wastewater,
    param_change_times,
    params::uciwweihr_model_params1;
    n_chains::Int64=1,
    seed::Int64=2024,
    verbose_optimize::Bool=false
)
    println("Getting init param values...")
    println("Using uciwweihr_model with wastewater.  Hardcoded sigma_ww and sigma_hosp!!!")
    obstimes_hosp = convert(Vector{Int64}, obstimes_hosp)
    obstimes_wastewater = convert(Vector{Int64}, obstimes_wastewater)
    param_change_times = convert(Vector{Int64}, param_change_times)

    my_model_optimize = uciwweihr_model(
        data_hosp,
        data_wastewater,
        obstimes_hosp,
        obstimes_wastewater,
        obstimes,
        param_change_times,
        params;
    )
    Random.seed!(seed)
    # Optimize
    MAP_init = optimize_many_MAP2(my_model_optimize, 10, 1, verbose_optimize)[1]
    Random.seed!(seed)
    MAP_noise = vcat(randn(length(MAP_init) - 1, n_chains), transpose(zeros(n_chains)))
    MAP_noise = [MAP_noise[:,i] for i in 1:size(MAP_noise,2)]
    init = repeat([MAP_init], n_chains) .+ 0.05 * MAP_noise
    return(init)
end



function uciwweihr_init_param(
    data_hosp,
    obstimes_hosp,
    param_change_times,
    params::uciwweihr_model_params3;
    n_chains::Int64=1,
    seed::Int64=2024,
    verbose_optimize::Bool=false
)
    println("Getting init param values...")
    println("Using uciwweihr_model w/ wastewater.  Hardcoded sigma_hosp and const hospitalization probability!!!")
    obstimes_hosp = convert(Vector{Int64}, obstimes_hosp)
    param_change_times = convert(Vector{Int64}, param_change_times)

    my_model_optimize = uciwweihr_model(
        data_hosp,
        obstimes_hosp,
        param_change_times,
        params;
    )
    Random.seed!(seed)
    # Optimize
    MAP_init = optimize_many_MAP2(my_model_optimize, 10, 1, verbose_optimize)[1]
    Random.seed!(seed)
    MAP_noise = vcat(randn(length(MAP_init) - 1, n_chains), transpose(zeros(n_chains)))
    MAP_noise = [MAP_noise[:,i] for i in 1:size(MAP_noise,2)]
    init = repeat([MAP_init], n_chains) .+ 0.05 * MAP_noise
    return(init)
end

