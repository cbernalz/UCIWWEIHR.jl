"""
    uciwweihr_sim_params
Struct for holding parameters used in the UCIWWEIHR ODE compartmental model simulation.

# Fields
- `time_points::Int64`: Number of time points for the simulation.
- `seed::Int64`: Seed for random number generation.
- `E_init::Int64`: Initial number of exposed individuals.
- `I_init::Int64`: Initial number of infected individuals.
- `H_init::Int64`: Initial number of hospitalized individuals.
- `gamma::Float64`: Rate of incubation.
- `nu::Float64`: Rate of leaving the infected compartment.
- `epsilon::Float64`: Rate of hospitalization recovery.
- `rho_gene::Float64`: Contribution of infected individual's pathogen genome into wastewater.
- `tau::Float64`: Scale/variation of the log concentration of pathogen genome in wastewater.
- `df::Float64`: Degrees of freedom for generalized t-distribution for log concentration of pathogen genome in wastewater.
- `sigma_hosp::Float64`: Standard deviation for the negative binomial distribution for hospital data.
- `Rt::Union{Float64, Vector{Float64}}`: Initial value or time series of the time-varying reproduction number.
- `sigma_Rt::Float64`: Standard deviation for random walk of time-varying reproduction number.
- `w::Union{Float64, Vector{Float64}}`: Initial value or time series of the time-varying hospitalization rate.
- `sigma_w::Float64`: Standard deviation for random walk of time-varying hospitalization rate.
"""
struct uciwweihr_sim_params
    time_points::Int64
    seed::Int64
    E_init::Int64
    I_init::Int64
    H_init::Int64
    gamma::Float64
    nu::Float64
    epsilon::Float64
    rho_gene::Float64
    tau::Float64
    df::Float64
    sigma_hosp::Float64
    Rt::Union{Float64, Vector{Float64}}
    sigma_Rt::Float64
    w::Union{Float64, Vector{Float64}}
    sigma_w::Float64
end

"""
    create_uciwweihr_params(; kwargs...)
Creates a `uciwweihr_sim_params` struct with the option to either use a predetermined `Rt` and `w` or generate them as random walks.

# Arguments
- `kwargs...`: Named arguments corresponding to the fields in `uciwweihr_sim_params`.

# Returns
- `params::uciwweihr_sim_params`: A struct with simulation parameters.
"""
function create_uciwweihr_params(; time_points::Int64=150, seed::Int64=1,
                                  E_init::Int64=200, I_init::Int64=100, H_init::Int64=20,
                                  gamma::Float64=1/4, nu::Float64=1/7, epsilon::Float64=1/5,
                                  rho_gene::Float64=0.011, tau::Float64=0.1, df::Float64=29.0,
                                  sigma_hosp::Float64=800.0,
                                  Rt::Union{Float64, Vector{Float64}}=1.0,
                                  sigma_Rt::Float64=sqrt(0.001),
                                  w::Union{Float64, Vector{Float64}}=0.35,
                                  sigma_w::Float64=sqrt(0.001))

    Random.seed!(seed)

    Rt_t = isa(Rt, Float64) ? generate_random_walk(time_points, sigma_Rt, Rt) : Rt
    w_t = isa(w, Float64) ? generate_random_walk(time_points, sigma_w, w) : w

    return uciwweihr_sim_params(time_points, seed, E_init, I_init, H_init, gamma, nu, epsilon,
                                rho_gene, tau, df, sigma_hosp, Rt_t, sigma_Rt, w_t, sigma_w)
end

"""
    generate_random_walk(time_points::Int64, sigma::Float64, init_val::Float64)
Generates a random walk time series.

# Arguments
- `time_points::Int64`: Number of time points.
- `sigma::Float64`: Standard deviation of the random walk.
- `init_val::Float64`: Initial value of the random walk.

# Returns
- `walk::Vector{Float64}`: Generated random walk.
"""
function generate_random_walk(time_points::Int64, sigma::Float64, init_val::Float64)
    walk = Float64[]
    log_val = log(init_val)
    for _ in 1:time_points
        log_val = rand(Normal(0, sigma)) + log_val
        push!(walk, exp(log_val))
    end
    return walk
end

"""
    generate_simulation_data(params::UCIWWEIHRParams)

Generates simulation data for the UCIWWEIHR ODE compartmental model.

# Arguments
- `params::uciwweihr_sim_params`: Struct containing parameters for the simulation.

# Returns
- `df::DataFrame`: A DataFrame containing the simulation data with columns `obstimes`, `log_ww_conc`, `hosp`, `rt`, and `wt`.
"""
function generate_simulation_data_uciwweihr(params::uciwweihr_sim_params)
    time_points = params.time_points

    alpha_t = params.Rt .* params.nu
    u0 = [params.E_init, params.I_init, params.H_init]
    p0 = [alpha_t[1], params.gamma, params.nu, params.w[1], params.epsilon]

    prob = ODEProblem(eihr_ode!, u0, (0.0, time_points), p0)
    
    function param_affect_beta!(integrator)
        ind_t = searchsortedfirst(1:time_points, integrator.t)
        integrator.p[1] = alpha_t[ind_t]
        integrator.p[4] = params.w[ind_t]
    end
    param_callback = PresetTimeCallback(1:time_points, param_affect_beta!, save_positions=(false, false))
    extra_ode_precision = false
    abstol = extra_ode_precision ? 1e-11 : 1e-9
    reltol = extra_ode_precision ? 1e-8 : 1e-6
    sol = solve(prob, Tsit5(); callback=param_callback, saveat=collect(1:time_points), save_start=true, 
                verbose=false, abstol=abstol, reltol=reltol, u0=u0, p=p0, tspan=(0.0, time_points))
    if sol.retcode != :Success
        Turing.@addlogprob! -Inf
        return
    end
    sol_array = Array(sol)
    I_comp_sol = clamp.(sol_array[2, 2:end], 1, 1e10)
    H_comp_sol = clamp.(sol_array[3, 2:end], 1, 1e10)

    # Log Gene Setup
    log_genes_mean = log.(I_comp_sol) .+ log(params.rho_gene)
    data_wastewater = [rand(GeneralizedTDist(log_genes_mean[t], params.tau, params.df)) for t in 1:time_points]
    data_hosp = [rand(NegativeBinomial2(H_comp_sol[t], params.sigma_hosp)) for t in 1:time_points]

    df = DataFrame(
        obstimes = 1:time_points,
        log_ww_conc = data_wastewater,
        hosp = data_hosp,
        rt = params.Rt,
        wt = params.w
    )

    return df
end
