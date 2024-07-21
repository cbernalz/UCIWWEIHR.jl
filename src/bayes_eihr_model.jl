# bayesian non-constant Rt and hosp rate model eihr
# -------------------------------------------------
"""
    bayes_eihr_model(...)
This is the bayesian semi-parametric model for the wastewater EIHR compartmental model.  
The defaults for this fuction will follow those of the default simulation in generate_simulation_data_ww_eihr.jl function.

# Arguments
- `data_hosp::Array{Float64}`: An array of hospital data.
- `data_wastewater::Array{Float64}`: An array of pathogen genome concentration in localized wastewater data.
- `obstimes::Array{Float64}`: An array of timepoints for observed hosp/wastewater.
- `E_init_sd::Float64`: Standard deviation for the initial number of exposed individuals.
- `E_init_mean::Int64`: Mean for the initial number of exposed individuals.
- `I_init_sd::Float64`: Standard deviation for the initial number of infected individuals.
- `I_init_mean::Int64`: Mean for the initial number of infected individuals.
- `H_init_sd::Float64`: Standard deviation for the initial number of hospitalized individuals.
- `H_init_mean::Int64`: Mean for the initial number of hospitalized individuals.
- `gamma_sd::Float64`: Standard deviation for the rate of incubation.
- `log_gamma_mean::Float64`: Mean for the rate of incubation on log scale.
- `nu_sd::Float64`: Standard deviation for the rate of laeving the infected compartment.
- `log_nu_mean::Float64`: Mean for the rate of leaving the infected compartment on the log scale.
- `epsilon_sd::Float64`: Standard deviation for the rate of hospitalization recovery.
- `log_epsilon_mean::Float64`: Mean for the rate of hospitalization recovery on the log scale.
- `rho_gene_sd::Float64`: Standard deviation for the rho prior.
- `log_rho_gene_mean::Float64`: Mean for the row prior on log scale.
- `tau_sd::Float64`: Standard deviation for the scale/variation of the log scale data.
- `log_tau_mean::Float64`: Mean for the scale/variation of the log scale data on log scale itself.
- `df_shape::Float64`: Shape parameter for the gamma distribution.
- `df_scale::Float64`: Scale parameter for the gamma distribution.
- `sigma_hosp_sd::Float64`: Standard deviation for the negative binomial distribution for hospital data.
- `sigma_hosp_mean::Float64`: Mean for the negative binomial distribution for hospital data.
- `Rt_init_sd::Float64`: Standard deviation for the initial value of the time-varying reproduction number.
- `Rt_init_mean::Float64`: Mean for the initial value of the time-varying reproduction number.
- `sigma_Rt_sd::Float64`: Standard deviation for normal prior of log time-varying reproduction number standard deviation.
- `sigma_Rt_mean::Float64`: Mean for normal prior of log time-varying reproduction number standard deviation.
- `w_init_sd::Float64`: Standard deviation for the initial value of the time-varying hospitalization rate.
- `w_init_mean::Float64`: Mean for the initial value of the time-varying hospitalization rate.
- `sigma_w_sd::Float64`: Standard deviation for normal prior of log time-varying hospitalization rate standard devaiation.
- `sigma_w_mean::Float64`: Mean for normal prior of time-varying hospitalization rate standard devaiation.
- `param_change_times::Array{Float64}`: An array of timepoints where the parameters change.

"""

@model function bayes_eihr_model(
    data_hosp::Array{Float64},
    data_wastewater::Array{Float64},
    obstimes::Array{Float64},
    E_init_sd::Float64, E_init_mean::Int64,
    I_init_sd::Float64, I_init_mean::Int64,
    H_init_sd::Float64, H_init_mean::Int64,
    gamma_sd::Float64, log_gamma_mean::Float64,
    nu_sd::Float64, log_nu_mean::Float64,
    epsilon_sd::Float64, log_epsilon_mean::Float64,
    rho_gene_sd::Float64, log_rho_gene_mean::Float64,
    tau_sd::Float64, log_tau_mean::Float64,
    df_shape::Float64, df_scale::Float64,
    sigma_hosp_sd::Float64, sigma_hosp_mean::Float64,
    Rt_init_sd::Float64, Rt_init_mean::Float64,
    sigma_Rt_sd::Float64, sigma_Rt_mean::Float64,
    w_init_sd::Float64, w_init_mean::Float64,
    sigma_w_sd::Float64, sigma_w_mean::Float64,
    param_change_times::Array{Float64}
    )


        # Prelims
        max_neg_bin_sigma = 1e10
        min_neg_bin_sigma = 1e-10


        # Calculate number of observed datapoints timepoints
        l_obs = length(obstimes)
        l_param_change_times = length(param_change_times)


        # PRIORS-----------------------------
        # Compartments
        E_init_non_centered ~ Normal()
        I_init_non_centered ~ Normal()
        H_init_non_centered ~ Normal()
        # Parameters for compartments
        gamma_non_centered ~ Normal()
        nu_non_centered ~ Normal()
        epsilon_non_centered ~ Normal()
        # Parameters for wastewater
        rho_gene_non_centered ~ Normal() # gene detection rate
        tau_non_centered ~ Normal() # standard deviation for log scale data
        df ~ Gamma(df_shape, df_scale)
        # Parameters for hospital
        sigma_hosp_non_centered ~ Normal()
        # Non-constant Rt
        Rt_params_non_centered ~ MvNormal(zeros(l_param_change_times + 2), I) # +2 for sigma and init
        sigma_Rt_non_centered = Rt_params_non_centered[1]
        Rt_init_non_centered = Rt_params_non_centered[2]
        log_Rt_steps_non_centered = Rt_params_non_centered[3:end]
        # Non-constant Hosp Rate w
        w_params_non_centered ~ MvNormal(zeros(l_param_change_times + 2), I) # +2 for sigma and init
        sigma_w_non_centered = w_params_non_centered[1]
        w_init_non_centered = w_params_non_centered[2]
        log_w_steps_non_centered = w_params_non_centered[3:end]


        # TRANSFORMATIONS--------------------
        # Compartments
        E_init = E_init_non_centered * E_init_sd + E_init_mean
        I_init = I_init_non_centered * I_init_sd + I_init_mean
        H_init = H_init_non_centered * H_init_sd + H_init_mean
        # Parameters for compartments
        gamma = exp(gamma_non_centered * gamma_sd + log_gamma_mean)
        nu = exp(nu_non_centered * nu_sd + log_nu_mean)
        epsilon = exp(epsilon_non_centered * epsilon_sd + log_epsilon_mean)
        # Parameters for wastewater
        rho_gene = exp(rho_gene_non_centered * rho_gene_sd + log_rho_gene_mean)
        tau = exp(tau_non_centered * tau_sd + log_tau_mean)
        # Parameters for hospital
        sigma_hosp = clamp.(sigma_hosp_non_centered * sigma_hosp_sd + sigma_hosp_mean, min_neg_bin_sigma, max_neg_bin_sigma)
        # Non-constant Rt
        Rt_init = exp(Rt_init_non_centered * Rt_init_sd + Rt_init_mean)
        sigma_Rt = exp(sigma_Rt_non_centered * sigma_Rt_sd + sigma_Rt_mean)
        alpha_t_no_init = exp.(log(Rt_init) .+ cumsum(log_Rt_steps_non_centered) * sigma_Rt) * nu
        alpha_init = Rt_init * nu
        alpha_t = vcat(alpha_init, alpha_t_no_init)
        # Non-constant Hosp Rate w
        w_init = exp(w_init_non_centered * w_init_sd + w_init_mean)
        sigma_w = exp(sigma_w_non_centered * sigma_w_sd + sigma_w_mean)
        w_no_init = exp.(log(w_init) .+ cumsum(log_w_steps_non_centered) * sigma_w)
        w_t = vcat(w_init, w_no_init)


        # ODE SETUP--------------------------
        prob = ODEProblem{true}(eihr_ode!, zeros(3), (0.0, obstimes[end]), ones(5))
        function param_affect_beta!(integrator)
            ind_t = searchsortedfirst(param_change_times, integrator.t) # Find the index of param_change_times that contains the current timestep
            integrator.p[1] = alpha_t_no_init[ind_t] # Replace alpha with a new value from alpha_t_no_init
            integrator.p[4] = w_no_init[ind_t] # Replace w with a new value from w_no_init
        end
        param_callback = PresetTimeCallback(param_change_times, param_affect_beta!, save_positions=(false, false))
        u0 = [E_init, I_init, H_init]
        p0 = [alpha_init, gamma, nu, w_init, epsilon]
        extra_ode_precision = false
        abstol = extra_ode_precision ? 1e-11 : 1e-9
        reltol = extra_ode_precision ? 1e-8 : 1e-6
        sol = solve(prob, Tsit5(); callback=param_callback, saveat=obstimes, save_start=true, 
                    verbose=false, abstol=abstol, reltol=reltol, u0=u0, p=p0, tspan=(0.0, obstimes[end]))
        # If the ODE solver fails, reject the sample by adding -Inf to the likelihood
        if sol.retcode != :Success
            println("An error occurred during ODE solution!!!")
            Turing.@addlogprob! -Inf
            return
        end
        sol_array = Array(sol)
        I_comp_sol = clamp.(sol_array[2,2:end],1, 1e10)


        # W-W means--------------------------
        # E - 1 // I - 2 // H - 3 // R - 4
        log_genes_mean = log.(I_comp_sol) .+ log(rho_gene) # first entry is the initial conditions, we want 2:end
        # Likelihood calculations------------
        sol_hosp = clamp.(sol_array[3,2:end], 1, 1e10)
        for i in 1:l_obs
            data_wastewater[i] ~ GeneralizedTDist(log_genes_mean[i], tau, df)
            data_hosp[i] ~ NegativeBinomial2(sol_hosp[i], sigma_hosp)
        end


        # Generated quantities
        H_comp = sol_array[3, :]
        rt_vals = alpha_t / nu
        w_t = w_t

        return (
            E_init,
            I_init,
            H_init,
            alpha_t = alpha_t,
            gamma = gamma,
            nu = nu,
            w_t = w_t,
            epsilon = epsilon,
            rt_vals = rt_vals,
            rho_gene = rho_gene,
            tau = tau,
            df = df,
            sigma_hosp = sigma_hosp,
            H = H_comp,
            log_genes_mean = log_genes_mean
        )


    end
