"""
    uciwweihr_likelihood_helpers(...)
    asdf

# Arguments
- asdf

"""
function uciwweihr_likelihood_helpers(
    obstimes_hosp,
    obstimes_wastewater,
    obstimes,
    param_change_times,
    params::uciwweihr_model_params2;
    E_init_non_centered, I_init_non_centered, H_init_non_centered,
    gamma_non_centered, nu_non_centered, epsilon_non_centered,
    rho_gene_non_centered, sigma_ww_non_centered, sigma_hosp_non_centered,
    Rt_params_non_centered, w_params_non_centered,
    warning_bool=true
)
    try
        # Prelims
        max_neg_bin_sigma = 1e10
        min_neg_bin_sigma = 1e-10

        # Non-constant Rt
        sigma_Rt_non_centered = Rt_params_non_centered[1]
        Rt_init_non_centered = Rt_params_non_centered[2]
        log_Rt_steps_non_centered = Rt_params_non_centered[3:end]
        # Non-constant Hosp Rate w
        sigma_w_non_centered = w_params_non_centered[1]
        w_init_non_centered = w_params_non_centered[2]
        logit_w_steps_non_centered = w_params_non_centered[3:end]


        # TRANSFORMATIONS--------------------
        # Compartments
        E_init = exp(E_init_non_centered * params.E_init_sd + params.log_E_init_mean)
        I_init = exp(I_init_non_centered * params.I_init_sd + params.log_I_init_mean)
        H_init = exp(H_init_non_centered * params.H_init_sd + params.log_H_init_mean)
        # Parameters for compartments
        gamma = exp(gamma_non_centered * params.gamma_sd + params.log_gamma_mean)
        nu = exp(nu_non_centered * params.nu_sd + params.log_nu_mean)
        epsilon = exp(epsilon_non_centered * params.epsilon_sd + params.log_epsilon_mean)
        # Parameters for wastewater
        rho_gene = exp(rho_gene_non_centered * params.rho_gene_sd + params.log_rho_gene_mean)
        sigma_ww = exp(sigma_ww_non_centered * params.sigma_ww_sd + params.log_sigma_ww_mean)
        # Parameters for hospital
        #sigma_hosp = clamp.(sigma_hosp_non_centered * params.sigma_hosp_sd + params.sigma_hosp_mean, min_neg_bin_sigma, max_neg_bin_sigma)    
        sigma_hosp = exp(sigma_hosp_non_centered * params.sigma_hosp_sd + params.sigma_hosp_mean)
        println("sigma_hosp = $sigma_hosp")
        println("pre_clamp = $(sigma_hosp_non_centered * params.sigma_hosp_sd + params.sigma_hosp_mean)")


        # Non-constant Rt
        Rt_init = exp(Rt_init_non_centered * params.Rt_init_sd + params.Rt_init_mean)
        sigma_Rt = exp(sigma_Rt_non_centered * params.sigma_Rt_sd + params.sigma_Rt_mean)
        alpha_t_no_init = exp.(log(Rt_init) .+ cumsum(log_Rt_steps_non_centered) * sigma_Rt) * nu
        alpha_init = Rt_init * nu
        alpha_t = vcat(alpha_init, alpha_t_no_init)
        Rt_t = alpha_t / nu
        # Non-constant Hosp Prob w
        w_init_logit = w_init_non_centered * params.w_init_sd + params.w_init_mean
        sigma_w = exp(sigma_w_non_centered * params.sigma_w_sd + params.sigma_w_mean)
        logit_w_no_init = w_init_logit .+ cumsum(logit_w_steps_non_centered) * sigma_w
        w_init = logistic(w_init_logit)
        w_no_init = logistic.(logit_w_no_init)
        w_t = vcat(w_init, w_no_init)

        # ODE SETUP--------------------------
        max_obstime_end = max(obstimes_hosp[end], obstimes_wastewater[end])
        prob = ODEProblem{true}(eihr_ode!, zeros(3), (0.0, max_obstime_end), ones(5))
        u0 = [E_init, I_init, H_init]
        p0 = (gamma, nu, epsilon, alpha_t, w_t, param_change_times)
        extra_ode_precision = false
        abstol = extra_ode_precision ? 1e-11 : 1e-9
        reltol = extra_ode_precision ? 1e-8 : 1e-6
        sol = solve(prob, Tsit5(); saveat=0.0:max_obstime_end, save_start=true, 
                    verbose=false, abstol=abstol, reltol=reltol, u0=u0, p=p0, tspan=(0.0, obstimes[end]))
        # If the ODE solver fails, reject the sample by adding -Inf to the likelihood
        if sol.retcode != :Success
            throw(ArgumentError("ODE solver failed!!!"))
        end

        sol_array = Array(sol)
        I_comp_sol = clamp.(sol_array[2,2:end],1, 1e10)
        E_comp_sol = clamp.(sol_array[1,2:end],1, 1e10)
        full_log_genes_mean = log.(I_comp_sol) .+ log(rho_gene) 
        H_comp_sol = clamp.(sol_array[3,2:end], 1, 1e10)
        H_means = H_comp_sol[obstimes_hosp]
        log_W_means = full_log_genes_mean[obstimes_wastewater]

        # Return --------------------------
        return (
            success = true,
            E_init = E_init, I_init = I_init, H_init = H_init,
            gamma = gamma, nu = nu, epsilon = epsilon,
            rho_gene = rho_gene, sigma_ww = sigma_ww, 
            sigma_hosp = sigma_hosp,
            Rt_t = Rt_t, Rt_init = Rt_init, sigma_Rt = sigma_Rt, alpha_t = alpha_t,
            w_init = w_init, sigma_w = sigma_w, w_t = w_t,
            E_comp_sol = E_comp_sol, I_comp_sol = I_comp_sol, H_comp_sol = H_comp_sol,
            H_means = H_means, log_W_means = log_W_means
        )
    catch e
        if warning_bool
            @warn "ODE solver or transformation failed: $e"
        end
        return (success = false,)
    end

end


function uciwweihr_likelihood_helpers(
    obstimes_hosp,
    obstimes_wastewater,
    obstimes,
    param_change_times,
    params::uciwweihr_model_params1;
    E_init_non_centered, I_init_non_centered, H_init_non_centered,
    gamma_non_centered, nu_non_centered, epsilon_non_centered,
    rho_gene_non_centered, 
    Rt_params_non_centered, w_params_non_centered,
    warning_bool=true
)
    try
        # Non-constant Rt
        sigma_Rt_non_centered = Rt_params_non_centered[1]
        Rt_init_non_centered = Rt_params_non_centered[2]
        log_Rt_steps_non_centered = Rt_params_non_centered[3:end]
        # Non-constant Hosp Rate w
        sigma_w_non_centered = w_params_non_centered[1]
        w_init_non_centered = w_params_non_centered[2]
        logit_w_steps_non_centered = w_params_non_centered[3:end]


        # TRANSFORMATIONS--------------------
        # Compartments
        E_init = exp(E_init_non_centered * params.E_init_sd + params.log_E_init_mean)
        I_init = exp(I_init_non_centered * params.I_init_sd + params.log_I_init_mean)
        H_init = exp(H_init_non_centered * params.H_init_sd + params.log_H_init_mean)
        # Parameters for compartments
        gamma = exp(gamma_non_centered * params.gamma_sd + params.log_gamma_mean)
        nu = exp(nu_non_centered * params.nu_sd + params.log_nu_mean)
        epsilon = exp(epsilon_non_centered * params.epsilon_sd + params.log_epsilon_mean)
        # Parameters for wastewater
        rho_gene = exp(rho_gene_non_centered * params.rho_gene_sd + params.log_rho_gene_mean)
        # Parameters for hospital 

        # Non-constant Rt
        Rt_init = exp(Rt_init_non_centered * params.Rt_init_sd + params.Rt_init_mean)
        sigma_Rt = exp(sigma_Rt_non_centered * params.sigma_Rt_sd + params.sigma_Rt_mean)
        alpha_t_no_init = exp.(log(Rt_init) .+ cumsum(log_Rt_steps_non_centered) * sigma_Rt) * nu
        alpha_init = Rt_init * nu
        alpha_t = vcat(alpha_init, alpha_t_no_init)
        Rt_t = alpha_t / nu
        # Non-constant Hosp Prob w
        w_init_logit = w_init_non_centered * params.w_init_sd + params.w_init_mean
        sigma_w = exp(sigma_w_non_centered * params.sigma_w_sd + params.sigma_w_mean)
        logit_w_no_init = w_init_logit .+ cumsum(logit_w_steps_non_centered) * sigma_w
        w_init = logistic(w_init_logit)
        w_no_init = logistic.(logit_w_no_init)
        w_t = vcat(w_init, w_no_init)


        # ODE SETUP--------------------------
        max_obstime_end = max(obstimes_hosp[end], obstimes_wastewater[end])
        prob = ODEProblem{true}(eihr_ode!, zeros(3), (0.0, max_obstime_end), ones(5))
        u0 = [E_init, I_init, H_init]
        p0 = (gamma, nu, epsilon, alpha_t, w_t, param_change_times)
        extra_ode_precision = false
        abstol = extra_ode_precision ? 1e-11 : 1e-9
        reltol = extra_ode_precision ? 1e-8 : 1e-6
        sol = solve(prob, Tsit5(); saveat=0.0:max_obstime_end, save_start=true, 
                    verbose=false, abstol=abstol, reltol=reltol, u0=u0, p=p0, tspan=(0.0, obstimes[end]))
        # If the ODE solver fails, reject the sample by adding -Inf to the likelihood
        if sol.retcode != :Success
            throw(ArgumentError("ODE solver failed!!!"))
        end
        sol_array = Array(sol)
        I_comp_sol = clamp.(sol_array[2,2:end],1, 1e10)
        E_comp_sol = clamp.(sol_array[1,2:end],1, 1e10)
        full_log_genes_mean = log.(I_comp_sol) .+ log(rho_gene) 
        H_comp_sol = clamp.(sol_array[3,2:end], 1, 1e10)
        H_means = H_comp_sol[obstimes_hosp]
        log_W_means = full_log_genes_mean[obstimes_wastewater]

        # Return --------------------------
        return (
            success = true,
            E_init = E_init, I_init = I_init, H_init = H_init,
            gamma = gamma, nu = nu, epsilon = epsilon,
            rho_gene = rho_gene, sigma_ww = params.sigma_wastewater, 
            sigma_hosp = params.sigma_hosp,
            Rt_t = Rt_t, Rt_init = Rt_init, sigma_Rt = sigma_Rt, alpha_t = alpha_t,
            w_init = w_init, sigma_w = sigma_w, w_t = w_t,
            E_comp_sol = E_comp_sol, I_comp_sol = I_comp_sol, H_comp_sol = H_comp_sol,
            H_means = H_means, log_W_means = log_W_means
        )
    catch e
        if warning_bool
            @warn "ODE solver or transformation failed: $e"
        end
        return (success = false,)
    end
end

function uciwweihr_likelihood_helpers(
    obstimes_hosp,
    param_change_times,
    params::uciwweihr_model_params2;
    E_init_non_centered, I_init_non_centered, H_init_non_centered,
    gamma_non_centered, nu_non_centered, epsilon_non_centered,
    sigma_hosp_non_centered,
    Rt_params_non_centered, w_params_non_centered,
    warning_bool=true
)
    try
        # Prelims
        max_neg_bin_sigma = 1e10
        min_neg_bin_sigma = 1e-10
        
        # Non-constant Rt
        sigma_Rt_non_centered = Rt_params_non_centered[1]
        Rt_init_non_centered = Rt_params_non_centered[2]
        log_Rt_steps_non_centered = Rt_params_non_centered[3:end]
        # Non-constant Hosp Rate w
        sigma_w_non_centered = w_params_non_centered[1]
        w_init_non_centered = w_params_non_centered[2]
        logit_w_steps_non_centered = w_params_non_centered[3:end]


        # TRANSFORMATIONS--------------------
        # Compartments
        E_init = exp(E_init_non_centered * params.E_init_sd + params.log_E_init_mean)
        I_init = exp(I_init_non_centered * params.I_init_sd + params.log_I_init_mean)
        H_init = exp(H_init_non_centered * params.H_init_sd + params.log_H_init_mean)
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
        Rt_t = alpha_t / nu
        # Non-constant Hosp Prob w
        w_init_logit = w_init_non_centered * params.w_init_sd + params.w_init_mean
        sigma_w = exp(sigma_w_non_centered * params.sigma_w_sd + params.sigma_w_mean)
        logit_w_no_init = w_init_logit .+ cumsum(logit_w_steps_non_centered) * sigma_w
        w_init = logistic(w_init_logit)
        w_no_init = logistic.(logit_w_no_init)
        w_t = vcat(w_init, w_no_init)


        # ODE SETUP--------------------------
        prob = ODEProblem{true}(eihr_ode!, zeros(3), (0.0, obstimes_hosp[end]), ones(5))
        u0 = [E_init, I_init, H_init]
        p0 = (gamma, nu, epsilon, alpha_t, w_t, param_change_times)
        extra_ode_precision = false
        abstol = extra_ode_precision ? 1e-11 : 1e-9
        reltol = extra_ode_precision ? 1e-8 : 1e-6
        sol = solve(prob, Tsit5(); saveat=0.0:obstimes_hosp[end], save_start=true, 
                    verbose=false, abstol=abstol, reltol=reltol, u0=u0, p=p0, tspan=(0.0, obstimes_hosp[end]))
        # If the ODE solver fails, reject the sample by adding -Inf to the likelihood
        if sol.retcode != :Success
            throw(ArgumentError("ODE solver failed!!!"))
        end
        sol_array = Array(sol)
        I_comp_sol = clamp.(sol_array[2,2:end],1, 1e10)
        E_comp_sol = clamp.(sol_array[1,2:end],1, 1e10)
        H_comp_sol = clamp.(sol_array[3,2:end], 1, 1e10)
        H_means = H_comp_sol[obstimes_hosp]

        # Return --------------------------
        return (
            success = true,
            E_init = E_init, I_init = I_init, H_init = H_init,
            gamma = gamma, nu = nu, epsilon = epsilon,
            sigma_hosp = sigma_hosp,
            Rt_t = Rt_t, Rt_init = Rt_init, sigma_Rt = sigma_Rt, alpha_t = alpha_t,
            w_init = w_init, sigma_w = sigma_w, w_t = w_t,
            E_comp_sol = E_comp_sol, I_comp_sol = I_comp_sol, H_comp_sol = H_comp_sol,
            H_means = H_means
        )
    catch e
        if warning_bool
            @warn "ODE solver or transformation failed: $e"
        end
        return (success = false,)
    end
end


function uciwweihr_likelihood_helpers(
    obstimes_hosp,
    param_change_times,
    params::uciwweihr_model_params1;
    E_init_non_centered, I_init_non_centered, H_init_non_centered,
    gamma_non_centered, nu_non_centered, epsilon_non_centered,
    Rt_params_non_centered, w_params_non_centered,
    warning_bool=true
)
    try
        # Non-constant Rt
        sigma_Rt_non_centered = Rt_params_non_centered[1]
        Rt_init_non_centered = Rt_params_non_centered[2]
        log_Rt_steps_non_centered = Rt_params_non_centered[3:end]
        # Non-constant Hosp Rate w
        sigma_w_non_centered = w_params_non_centered[1]
        w_init_non_centered = w_params_non_centered[2]
        logit_w_steps_non_centered = w_params_non_centered[3:end]


        # TRANSFORMATIONS--------------------
        # Compartments
        E_init = exp(E_init_non_centered * params.E_init_sd + params.log_E_init_mean)
        I_init = exp(I_init_non_centered * params.I_init_sd + params.log_I_init_mean)
        H_init = exp(H_init_non_centered * params.H_init_sd + params.log_H_init_mean)
        # Parameters for compartments
        gamma = exp(gamma_non_centered * params.gamma_sd + params.log_gamma_mean)
        nu = exp(nu_non_centered * params.nu_sd + params.log_nu_mean)
        epsilon = exp(epsilon_non_centered * params.epsilon_sd + params.log_epsilon_mean)

        # Non-constant Rt
        Rt_init = exp(Rt_init_non_centered * params.Rt_init_sd + params.Rt_init_mean)
        sigma_Rt = exp(sigma_Rt_non_centered * params.sigma_Rt_sd + params.sigma_Rt_mean)
        alpha_t_no_init = exp.(log(Rt_init) .+ cumsum(log_Rt_steps_non_centered) * sigma_Rt) * nu
        alpha_init = Rt_init * nu
        alpha_t = vcat(alpha_init, alpha_t_no_init)
        Rt_t = alpha_t / nu
        # Non-constant Hosp Prob w
        w_init_logit = w_init_non_centered * params.w_init_sd + params.w_init_mean
        sigma_w = exp(sigma_w_non_centered * params.sigma_w_sd + params.sigma_w_mean)
        logit_w_no_init = w_init_logit .+ cumsum(logit_w_steps_non_centered) * sigma_w
        w_init = logistic(w_init_logit)
        w_no_init = logistic.(logit_w_no_init)
        w_t = vcat(w_init, w_no_init)


        # ODE SETUP--------------------------
        prob = ODEProblem{true}(eihr_ode!, zeros(3), (0.0, obstimes_hosp[end]), ones(5))
        u0 = [E_init, I_init, H_init]
        p0 = (gamma, nu, epsilon, alpha_t, w_t, param_change_times)
        extra_ode_precision = false
        abstol = extra_ode_precision ? 1e-11 : 1e-9
        reltol = extra_ode_precision ? 1e-8 : 1e-6
        sol = solve(prob, Tsit5(); saveat=0.0:obstimes_hosp[end], save_start=true, 
                    verbose=false, abstol=abstol, reltol=reltol, u0=u0, p=p0, tspan=(0.0, obstimes_hosp[end]))
        # If the ODE solver fails, reject the sample by adding -Inf to the likelihood
        if sol.retcode != :Success
            throw(ArgumentError("ODE solver failed!!!"))
        end
        sol_array = Array(sol)
        I_comp_sol = clamp.(sol_array[2,2:end],1, 1e10)
        E_comp_sol = clamp.(sol_array[1,2:end],1, 1e10)
        H_comp_sol = clamp.(sol_array[3,2:end], 1, 1e10)
        H_means = H_comp_sol[obstimes_hosp]

        # Return --------------------------
        return (
            success = true,
            E_init = E_init, I_init = I_init, H_init = H_init,
            gamma = gamma, nu = nu, epsilon = epsilon,
            sigma_hosp = params.sigma_hosp,
            Rt_t = Rt_t, Rt_init = Rt_init, sigma_Rt = sigma_Rt, alpha_t = alpha_t,
            w_init = w_init, sigma_w = sigma_w, w_t = w_t,
            E_comp_sol = E_comp_sol, I_comp_sol = I_comp_sol, H_comp_sol = H_comp_sol,
            H_means = H_means
        )
    catch e
        if warning_bool
            @warn "ODE solver or transformation failed: $e"
        end
        return (success = false,)
    end
end