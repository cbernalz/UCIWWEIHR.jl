"""
## Generating Simulation Data for Wastewater EIHR

To generate simulation data for wastewater EIHR, you can use the `generate_simulation_data_ww_eihr` function defined in the `UCIWWEIHR.jl` package. This function allows you to customize various parameters for the simulation.

### Function Signature

```julia

# Arguments
- time_points::Int64: Number of time points wanted for simulation. Default value is 150.
- seed::Int64: Seed for random number generation. Default value is 1.
- E_init::Int64: Initial number of exposed individuals. Default value is 200.
- I_init::Int64: Initial number of infected individuals. Default value is 100.
- H_init::Int64: Initial number of hospitalized individuals. Default value is 20.
- gamma::Float64: Rate of incubation. Default value is 1/4.
- nu::Float64: Rate of leaving the infected compartment. Default value is 1/7.
- epsilon::Float64: Rate of hospitalization recovery. Default value is 1/5.
- rho_gene::Float64: Contribution of infected individual's pathogen genome into wastewater. Default value is 0.011.
- tau::Float64: Scale/variation of the log concentration of pathogen genome in wastewater. Default value is 0.1.
- df::Float64: Degrees of freedom for generalized t distribution for log concentration of pathogen genome in wastewater. Default value is 29.
- sigma_hosp::Float64: Standard deviation for the negative binomial distribution for hospital data. Default value is 800.
- Rt_init::Float64: Initial value of the time-varying reproduction number. Default value is 1.
- sigma_Rt::Float64: Standard deviation for random walk of time-varying reproduction number. Default value is sqrt(0.02).
- w_init::Float64: Initial value of the time-varying hospitalization rate. Default value is 0.35.
- sigma_w::Float64: Standard deviation for random walk of time-varying hospitalization rate. Default value is sqrt(0.02).
"""
function generate_simulation_data_ww_eihr(
    time_points::Int64=150, seed::Int64=1,
    E_init::Int64=200, I_init::Int64=100, H_init::Int64=20,
    gamma::Float64=1/4, nu::Float64=1/7, epsilon::Float64=1/5,
    rho_gene::Float64=0.011, tau::Float64=0.1, df::Float64=29.0,
    sigma_hosp::Float64=800.0,
    Rt_init::Float64=1.0, sigma_Rt::Float64=sqrt(0.001),
    w_init::Float64=0.35, sigma_w::Float64=sqrt(0.001),
)
    
    Random.seed!(seed)
    
    # Rt and W SETUP--------------------------
    Rt_t_no_init = Float64[]  # Pre-defined vector
    w_no_init = Float64[]  # Pre-defined vector
    log_Rt_t = log(Rt_init)
    log_w_t = log(w_init)
    for _ in 1:time_points
        log_Rt_t = log_Rt_t + rand(Normal(0, sigma_Rt))
        log_w_t = log_w_t + rand(Normal(0, sigma_w))
        push!(Rt_t_no_init, exp(log_Rt_t))
        push!(w_no_init, exp(log_w_t))
    end
    alpha_t_no_init = Rt_t_no_init * nu
    alpha_init = Rt_init * nu

    # ODE SETUP--------------------------
    prob = ODEProblem{true}(eihr_ode!, zeros(3), (0.0, time_points), ones(5))
    function param_affect_beta!(integrator)
        ind_t = searchsortedfirst(collect(1:time_points), integrator.t) # Find the index of collect(1:time_points) that contains the current timestep
        integrator.p[1] = alpha_t_no_init[ind_t] # Replace alpha with a new value from alpha_t_no_init
        integrator.p[4] = w_no_init[ind_t] # Replace w with a new value from w_no_init
    end
    param_callback = PresetTimeCallback(collect(1:time_points), param_affect_beta!, save_positions=(false, false))
    u0 = [E_init, I_init, H_init]
    p0 = [alpha_init, gamma, nu, w_init, epsilon]
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
    I_comp_sol = clamp.(sol_array[2,2:end],1, 1e10)
    H_comp_sol = clamp.(sol_array[3,2:end], 1, 1e10)

    # Log Gene SETUP--------------------------
    log_genes_mean = log.(I_comp_sol) .+ log(rho_gene) # first entry is the initial conditions, we want 2:end
    data_wastewater = zeros(time_points)
    data_hosp = zeros(time_points)
    for t_i in 1:time_points
        data_wastewater[t_i] = rand(GeneralizedTDist(log_genes_mean[t_i], tau, df))
        data_hosp[t_i] = rand(NegativeBinomial2(H_comp_sol[t_i], sigma_hosp))
    end

    df = DataFrame(
        obstimes = 1:time_points,
        log_ww_conc = data_wastewater,
        hosp = data_hosp,
        rt = Rt_t_no_init
    );
    return df

    
end
