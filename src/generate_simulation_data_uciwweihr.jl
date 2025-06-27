"""
    uciwweihr_sim_params

Struct for holding parameters used in the UCIWWEIHR ODE compartmental model simulation.  Use `create_uciwweihr_sim_params` to create an instance of this struct.

# Fields
- `time_points::Int64`: Number of time points for the simulation.
- `seed::Int64`: Seed for random number generation.
- `E_init::Int64`: Initial number of exposed individuals.
- `I_init::Int64`: Initial number of infected individuals.
- `H_init::Int64`: Initial number of hospitalized individuals.
- `H_init_inc::Int64`: Initial number of hospitalized individuals for incidence.
- `gamma::Float64`: Rate of incubation.
- `nu::Float64`: Rate of leaving the infected compartment.
- `epsilon::Float64`: Rate of hospitalization recovery.
- `rho_gene::Float64`: Contribution of infected individual's pathogen genome into wastewater.
- `sigma_ww::Float64`: standard deviation of the log concentration of pathogen genome in wastewater.
- `sigma_hosp::Float64`: Overdispersion for the negative binomial distribution for hospital data.
- `sigma_hosp_inc::Float64`: Overdispersion for the negative binomial distribution for hospital incidence data.
- `Rt::Union{Float64, Vector{Float64}}`: Initial value or time series of the time-varying reproduction number.
- `sigma_Rt::Float64`: Standard deviation for random walk of time-varying reproduction number.
- `w::Union{Float64, Vector{Float64}}`: Initial value or time series of the time-varying hospitalization rate.
- `sigma_w::Float64`: Standard deviation for random walk of time-varying hospitalization rate.
- `rt_init::Float64`: Initial value of the time-varying reproduction number, NOT USER SPECIFIED `create_uciwweihr_params` TAKES CARE OF THIS.
- `w_init::Float64`: Initial value of the time-varying hospitalization rate, NOT USER SPECIFIED `create_uciwweihr_params` TAKES CARE OF THIS.
"""
struct uciwweihr_sim_params
    time_points::Union{Int64, Nothing}
    seed::Union{Int64, Nothing}
    E_init::Union{Int64, Nothing}
    I_init::Union{Int64, Nothing}
    H_init::Union{Int64, Nothing}
    H_init_inc::Union{Int64, Nothing}
    gamma::Union{Float64, Nothing}
    nu::Union{Float64, Nothing}
    epsilon::Union{Float64, Nothing}
    rho_gene::Union{Float64, Nothing}
    sigma_ww::Union{Float64, Nothing}
    sigma_hosp::Union{Float64, Nothing}
    sigma_hosp_inc::Union{Float64, Nothing}
    Rt::Union{Float64, Vector{Float64}, Nothing}
    sigma_Rt::Union{Float64, Nothing}
    w::Union{Float64, Vector{Float64}, Nothing}
    sigma_w::Union{Float64, Nothing}
    rt_init::Union{Float64, Nothing}
    w_init::Union{Float64, Nothing}
end


"""
create_uciwweihr_sim_params(; kwargs...)

Creates a `uciwweihr_sim_params` struct with the option to either use a predetermined `Rt` and `w` or generate them as random walks.

# Arguments
- `kwargs...`: Named arguments corresponding to the fields in `uciwweihr_sim_params`.

# Returns
- `params::uciwweihr_sim_params`: A struct with simulation parameters.
"""
function create_uciwweihr_sim_params(; time_points::Int64=150, seed::Int64=1,
                                  E_init::Int64=200, I_init::Int64=100, H_init::Int64=20, H_init_inc::Int64=5,
                                  gamma::Float64=1/4, nu::Float64=1/7, epsilon::Float64=1/5,
                                  rho_gene::Float64=0.011, sigma_ww::Float64=0.2,
                                  sigma_hosp::Float64=350.0,
                                  sigma_hosp_inc::Float64=350.0,
                                  Rt::Union{Float64, Vector{Float64}}=1.0,
                                  sigma_Rt::Float64=sqrt(0.001),
                                  w::Union{Float64, Vector{Float64}}=0.35,
                                  sigma_w::Float64=sqrt(0.001))

    Random.seed!(seed)
    rt_init = isa(Rt, Float64) ? Rt : Rt[1]
    w_init = isa(w, Float64) ? w : w[1]

    Rt_t = isa(Rt, Float64) ? generate_random_walk(time_points, sigma_Rt, Rt, seed) : Rt
    w_t = isa(w, Float64) ? generate_logit_normal_random_walk(time_points, sigma_w, w, seed) : w

    return uciwweihr_sim_params(time_points, seed, E_init, I_init, H_init, H_init_inc, gamma, nu, epsilon,
                                rho_gene, sigma_ww, sigma_hosp, sigma_hosp_inc, Rt_t, sigma_Rt, w_t, sigma_w, rt_init, w_init)
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
function generate_random_walk(time_points::Int64, sigma::Float64, init_val::Float64, seed)
    Random.seed!(seed)
    walk = Float64[]
    log_val = log(init_val)
    for _ in 1:time_points
        log_val = rand(Normal(0, sigma)) + log_val
        push!(walk, exp(log_val))
    end
    return walk
end

"""
    generate_logit_normal_random_walk(time_points::Int64, sigma::Float64, init_val::Float64)

Generates a logit-normal random walk time series.

# Arguments
- `time_points::Int64`: Number of time points.
- `sigma::Float64`: Standard deviation of the random walk in logit space.
- `init_val::Float64`: Initial value of the random walk on the probability scale.

# Returns
- `walk::Vector{Float64}`: Generated random walk on the probability scale.
"""
function generate_logit_normal_random_walk(time_points::Int64, sigma::Float64, init_val::Float64, seed)
    Random.seed!(seed)
    walk = Float64[]
    logit_val = logit(init_val)
    for _ in 1:time_points
        logit_val = rand(Normal(0, sigma)) + logit_val
        push!(walk, logistic(logit_val))
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
    Random.seed!(params.seed)
    time_points = params.time_points

    # Preping weekly params
    alpha_t = params.Rt .* params.nu

    # ODE SETUP--------------------------
    u0 = [params.E_init, params.I_init, params.H_init, params.H_init_inc]
    p0 = (params.gamma, params.nu, params.epsilon, alpha_t, params.w, 1:time_points)
    prob = ODEProblem(eihr_ode!, u0, (1.0, time_points), p0)
    extra_ode_precision = false
    abstol = extra_ode_precision ? 1e-11 : 1e-9
    reltol = extra_ode_precision ? 1e-8 : 1e-6
    sol = solve(prob, Tsit5(); saveat=1:time_points, save_start=true, 
                verbose=false, abstol=abstol, reltol=reltol, u0=u0, p=p0, tspan=(1, time_points))
    # If the ODE solver fails, reject the sample by adding -Inf to the likelihood
    if sol.retcode != :Success
        Turing.@addlogprob! -Inf
        return
    end
    sol_array = Array(sol)
    I_comp_sol = clamp.(sol_array[2, 1:end], 1, 1e10)
    H_comp_sol = clamp.(sol_array[3, 1:end], 1, 1e10)
    pushfirst!(sol_array[4, 2:end] - sol_array[4, 1:end-1], params.H_init_inc)
    H_inc_comp_sol = clamp.(pushfirst!(sol_array[4, 2:end] - sol_array[4, 1:end-1], params.H_init_inc), 1, 1e10)

    # Log Gene Setup
    log_genes_mean = log.(I_comp_sol) .+ log(params.rho_gene)

    # Data
    data_wastewater = [rand(Normal(log_genes_mean[t], params.sigma_ww)) for t in 1:time_points]
    data_hosp = [rand(NegativeBinomial2(H_comp_sol[t], params.sigma_hosp)) for t in 1:time_points]
    data_hosp_inc = [rand(NegativeBinomial2(H_inc_comp_sol[t], params.sigma_hosp_inc)) for t in 1:time_points]

    df = DataFrame(
        obstimes = 1:time_points,
        log_ww_conc = data_wastewater,
        hosp = data_hosp,
        hosp_inc = data_hosp_inc,
        rt = params.Rt,
        wt = params.w,
        E_ode_comp_sol = clamp.(sol_array[1, 1:end], 1, 1e10),
        I_ode_comp_sol = I_comp_sol,
        H_ode_comp_sol = H_comp_sol,
        H_ode_inc_comp_sol = H_inc_comp_sol,
    )

    return df
end
