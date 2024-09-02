"""
    uciwweihr_model(...)
This is the bayesian semi-parametric model for the wastewater EIHR compartmental model.  
The defaults for this fuction will follow those of the default simulation in generate_simulation_data_ww_eihr.jl function.

# Arguments
- `data_hosp`: An array of hospital data.
- `data_wastewater`: An array of pathogen genome concentration in localized wastewater data.  If this is not avaliable, the model used will be one that only uses hospital data.
- `obstimes_hosp`: An array of timepoints for observed hospital data.
- `obstimes_wastewater`: An array of timepoints for observed wastewater data.
- `param_change_times`: An array of timepoints where the parameters change.
- `params::uciwweihr_model_params`: A struct containing parameters for the model.

"""
@model function uciwweihr_model(
    data_hosp,
    data_wastewater,
    obstimes_hosp,
    obstimes_wastewater;
    param_change_times,
    params::uciwweihr_model_params
    )

        # Prelims
        max_neg_bin_sigma = 1e10
        min_neg_bin_sigma = 1e-10


        # Calculate number of observed datapoints timepoints
        l_obs_hosp = length(obstimes_hosp)
        l_obs_ww = length(obstimes_wastewater)
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
        df ~ Gamma(params.df_shape, params.df_scale)
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
        logit_w_steps_non_centered = w_params_non_centered[3:end]


        # TRANSFORMATIONS--------------------
        # Compartments
        E_init = E_init_non_centered * params.E_init_sd + params.E_init_mean
        I_init = I_init_non_centered * params.I_init_sd + params.I_init_mean
        H_init = H_init_non_centered * params.H_init_sd + params.H_init_mean
        # Parameters for compartments
        gamma = exp(gamma_non_centered * params.gamma_sd + params.log_gamma_mean)
        nu = exp(nu_non_centered * params.nu_sd + params.log_nu_mean)
        epsilon = exp(epsilon_non_centered * params.epsilon_sd + params.log_epsilon_mean)
        # Parameters for wastewater
        rho_gene = exp(rho_gene_non_centered * params.rho_gene_sd + params.log_rho_gene_mean)
        tau = exp(tau_non_centered * params.tau_sd + params.log_tau_mean)
        # Parameters for hospital
        sigma_hosp = clamp.(sigma_hosp_non_centered * params.sigma_hosp_sd + params.sigma_hosp_mean, min_neg_bin_sigma, max_neg_bin_sigma)
        # Non-constant Rt
        Rt_init = exp(Rt_init_non_centered * params.Rt_init_sd + params.Rt_init_mean)
        sigma_Rt = exp(sigma_Rt_non_centered * params.sigma_Rt_sd + params.sigma_Rt_mean)
        alpha_t_no_init = exp.(log(Rt_init) .+ cumsum(log_Rt_steps_non_centered) * sigma_Rt) * nu
        alpha_init = Rt_init * nu
        alpha_t = vcat(alpha_init, alpha_t_no_init)
        # Non-constant Hosp Prob w
        w_init_logit = w_init_non_centered * params.w_init_sd + params.w_init_mean
        sigma_w = exp(sigma_w_non_centered * params.sigma_w_sd + params.sigma_w_mean)
        logit_w_no_init = w_init_logit .+ cumsum(logit_w_steps_non_centered) * sigma_w
        w_init = logistic(w_init_logit)
        w_no_init = logistic.(logit_w_no_init)
        w_t = vcat(w_init, w_no_init)


        # ODE SETUP--------------------------
        max_obstime_end = max(obstimes_hosp[end], obstimes_wastewater[end])
        obstimes = unique(vcat(obstimes_hosp, obstimes_wastewater))
        obstimes = sort(obstimes)
        prob = ODEProblem{true}(eihr_ode!, zeros(3), (0.0, max_obstime_end), ones(5))
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
        sol = solve(prob, Tsit5(); callback=param_callback, saveat=0.0:max_obstime_end, save_start=true, 
                    verbose=false, abstol=abstol, reltol=reltol, u0=u0, p=p0, tspan=(0.0, obstimes[end]))
        # If the ODE solver fails, reject the sample by adding -Inf to the likelihood
        if sol.retcode != :Success
            Turing.@addlogprob! -Inf
            return
        end
        obstimes_hosp_indices = Int.(obstimes_hosp)
        obstimes_wastewater_indices = Int.(obstimes_wastewater)
        sol_array = Array(sol)
        I_comp_sol = clamp.(sol_array[2,2:end],1, 1e10)
        full_log_genes_mean = log.(I_comp_sol) .+ log(rho_gene) 
        H_comp_sol = clamp.(sol_array[3,2:end], 1, 1e10)
        H_means = H_comp_sol[obstimes_hosp_indices]
        log_W_means = full_log_genes_mean[obstimes_wastewater_indices]


        # W-W means--------------------------
        # E - 1 // I - 2 // H - 3 // R - 4
        # Likelihood calculations------------
        for i in 1:l_obs_ww
            data_wastewater[i] ~ GeneralizedTDist(log_W_means[i], tau, df)
        end
        for i in 1:l_obs_hosp
            data_hosp[i] ~ NegativeBinomial2(H_means[i], sigma_hosp)
        end


        

        # Generated quantities
        rt_vals = alpha_t_no_init / nu
        rt_init = alpha_init / nu
        w_t = w_no_init

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
            H = H_means,
            log_genes_mean = log_W_means,
            rt_init = rt_init,
            w_init = w_init
        )


    end



@model function uciwweihr_model(
    data_hosp,
    obstimes_hosp;
    param_change_times,
    params::uciwweihr_model_params
    )
    
    
        # Prelims
        max_neg_bin_sigma = 1e10
        min_neg_bin_sigma = 1e-10
    
    
        # Calculate number of observed datapoints timepoints
        l_obs = length(obstimes_hosp)
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
        logit_w_steps_non_centered = w_params_non_centered[3:end]
    
    
        # TRANSFORMATIONS--------------------
        # Compartments
        E_init = E_init_non_centered * params.E_init_sd + params.E_init_mean
        I_init = I_init_non_centered * params.I_init_sd + params.I_init_mean
        H_init = H_init_non_centered * params.H_init_sd + params.H_init_mean
        # Parameters for compartments
        gamma = exp(gamma_non_centered * params.gamma_sd + params.log_gamma_mean)
        nu = exp(nu_non_centered * params.nu_sd + params.log_nu_mean)
        epsilon = exp(epsilon_non_centered * params.epsilon_sd + params.log_epsilon_mean)
        # Parameters for hospital
        sigma_hosp = clamp.(sigma_hosp_non_centered * params.sigma_hosp_sd + params.sigma_hosp_mean, min_neg_bin_sigma, max_neg_bin_sigma)
        # Non-constant Rt
        Rt_init = exp(Rt_init_non_centered * params.Rt_init_sd + params.Rt_init_mean)
        sigma_Rt = exp(sigma_Rt_non_centered * params.sigma_Rt_sd + params.sigma_Rt_mean)
        alpha_t_no_init = exp.(log(Rt_init) .+ cumsum(log_Rt_steps_non_centered) * sigma_Rt) * nu
        alpha_init = Rt_init * nu
        alpha_t = vcat(alpha_init, alpha_t_no_init)
        # Non-constant Hosp Prob w
        w_init_logit = w_init_non_centered * params.w_init_sd + params.w_init_mean
        sigma_w = exp(sigma_w_non_centered * params.sigma_w_sd + params.sigma_w_mean)
        logit_w_no_init = w_init_logit .+ cumsum(logit_w_steps_non_centered) * sigma_w
        w_init = logistic(w_init_logit)
        w_no_init = logistic.(logit_w_no_init)
        w_t = vcat(w_init, w_no_init)
    
    
        # ODE SETUP--------------------------
        prob = ODEProblem{true}(eihr_ode!, zeros(3), (0.0, obstimes_hosp[end]), ones(5))
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
        sol = solve(prob, Tsit5(); callback=param_callback, saveat=obstimes_hosp, save_start=true, 
                    verbose=false, abstol=abstol, reltol=reltol, u0=u0, p=p0, tspan=(0.0, obstimes_hosp[end]))
        # If the ODE solver fails, reject the sample by adding -Inf to the likelihood
        if sol.retcode != :Success
            Turing.@addlogprob! -Inf
            return
        end
        sol_array = Array(sol)
        I_comp_sol = clamp.(sol_array[2,2:end],1, 1e10)
    
    
        # Likelihood calculations------------
        sol_hosp = clamp.(sol_array[3,2:end], 1, 1e10)
        for i in 1:l_obs
            data_hosp[i] ~ NegativeBinomial2(sol_hosp[i], sigma_hosp)
        end
    
    
        # Generated quantities
        rt_vals = alpha_t_no_init / nu
        rt_init = alpha_init / nu
        w_t = w_no_init
    
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
            sigma_hosp = sigma_hosp,
            H = H_comp,
            rt_init = rt_init,
            w_init = w_init
        )
    
    
    end
    