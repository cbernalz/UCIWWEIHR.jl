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
- `return_bool`: A boolean to indicate if the model is to use the return statement.  **Only set to false if only forecast is desired**

"""
@model function uciwweihr_model(
    data_hosp,
    data_wastewater,
    obstimes_hosp,
    obstimes_wastewater,
    param_change_times,
    params::model_params_time_var_hosp;
    incidence_model_bool=false,
    warning_bool=true
    )
        # hosp_ww model w/ time-varying w

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
        rho_gene_non_centered ~ Normal()
        sigma_ww_non_centered ~ Normal()
        # Parameters for hospital
        sigma_hosp_non_centered ~ Normal()
        # Non-constant Rt
        Rt_params_non_centered ~ MvNormal(zeros(length(param_change_times) - 1 + 2), I) # +2 for sigma and init
        # Non-constant Hosp Rate w
        w_params_non_centered ~ MvNormal(zeros(length(param_change_times) - 1 + 2), I) # +2 for sigma and init

        # TRANSFORMATIONS-----------------------------
        trans = likelihood_helpers(
            obstimes_hosp,
            obstimes_wastewater,
            param_change_times,
            params;
            incidence_model=incidence_model_bool,
            E_init_non_centered, I_init_non_centered, H_init_non_centered,
            gamma_non_centered, nu_non_centered, epsilon_non_centered,
            rho_gene_non_centered, sigma_ww_non_centered, sigma_hosp_non_centered,
            Rt_params_non_centered, w_params_non_centered,
            warning_bool=warning_bool
        )
        # Reject if the helper function failed and skip sample
        if !trans.success
            Turing.@addlogprob! -Inf
            return
        end

        # Likelihood calculations------------
        for i in 1:length(obstimes_wastewater)
            data_wastewater[i] ~ Normal(trans.log_W_means[i], trans.sigma_ww)
        end
        for i in 1:length(obstimes_hosp)
            if incidence_model_bool
                # If incidence model, use the incidence data
                data_hosp[i] ~ NegativeBinomial2(trans.H_inc_means[i], trans.sigma_hosp)
            else
                data_hosp[i] ~ NegativeBinomial2(trans.H_prev_means[i], trans.sigma_hosp)
            end
        end

        if incidence_model_bool
            return (
                E_init = trans.E_init, I_init = trans.I_init, H_init = trans.H_init,
                alpha_t = trans.alpha_t, w_t = trans.w_t, rt_vals = trans.Rt_t,
                gamma = trans.gamma, nu = trans.nu, epsilon = trans.epsilon,
                sigma_Rt = trans.sigma_Rt, sigma_w = trans.sigma_w, rt_init = trans.Rt_init, w_init = trans.w_init,
                rho_gene = trans.rho_gene,
                sigma_ww = trans.sigma_ww, sigma_hosp = trans.sigma_hosp,
                H = trans.H_comp_sol, I = trans.I_comp_sol, E = trans.E_comp_sol, CH = trans.CH_comp_sol, H_inc_comp_sol = trans.H_inc_comp_sol,
                H_inc_means = trans.H_inc_means, log_genes_mean = trans.log_W_means,
            )
        else
            return (
                E_init = trans.E_init, I_init = trans.I_init, H_init = trans.H_init,
                alpha_t = trans.alpha_t, w_t = trans.w_t, rt_vals = trans.Rt_t,
                gamma = trans.gamma, nu = trans.nu, epsilon = trans.epsilon,
                sigma_Rt = trans.sigma_Rt, sigma_w = trans.sigma_w, rt_init = trans.Rt_init, w_init = trans.w_init,
                rho_gene = trans.rho_gene,
                sigma_ww = trans.sigma_ww, sigma_hosp = trans.sigma_hosp,
                H = trans.H_comp_sol, I = trans.I_comp_sol, E = trans.E_comp_sol, CH = trans.CH_comp_sol, H_inc_comp_sol = trans.H_inc_comp_sol,
                H_prev_means = trans.H_prev_means, log_genes_mean = trans.log_W_means,
            )
        end

    end


@model function uciwweihr_model(
    data_hosp,
    obstimes_hosp,
    param_change_times,
    params::model_params_non_time_var_hosp_no_ww;
    warning_bool=true
    )
        # hosp_only model without time-varying hosp probability

        # PRIORS-----------------------------
        # Compartments
        E_init_non_centered ~ Normal()
        I_init_non_centered ~ Normal()
        H_init_non_centered ~ Normal()
        # Parameters for compartments
        gamma_non_centered ~ Normal()
        nu_non_centered ~ Normal()
        epsilon_non_centered ~ Normal()
        w_param_non_centered ~ Normal()
        # Parameters for hospital
        sigma_hosp_non_centered ~ Normal()
        # Non-constant Rt
        Rt_params_non_centered ~ MvNormal(zeros(length(param_change_times) - 1 + 2), I) # +2 for sigma and init

        # TRANSFORMATIONS-----------------------------
        trans = likelihood_helpers(
            obstimes_hosp,
            param_change_times,
            params;
            E_init_non_centered, I_init_non_centered, H_init_non_centered,
            gamma_non_centered, nu_non_centered, epsilon_non_centered,
            sigma_hosp_non_centered,
            Rt_params_non_centered, w_param_non_centered,
            warning_bool=warning_bool
        )
        # Reject if the helper function failed and skip sample
        if !trans.success
            Turing.@addlogprob! -Inf
            return
        end

        # Likelihood calculations------------
        for i in 1:length(obstimes_hosp)
            data_hosp[i] ~ NegativeBinomial2(trans.H_prev_means[i], trans.sigma_hosp)
        end

        return (
            E_init = trans.E_init,
            I_init = trans.I_init,
            H_init = trans.H_init,
            alpha_t = trans.alpha_t,
            gamma = trans.gamma,
            nu = trans.nu,
            epsilon = trans.epsilon,
            w = trans.w,
            rt_vals = trans.Rt_t,
            sigma_Rt = trans.sigma_Rt,
            sigma_hosp = trans.sigma_hosp,
            H = trans.H_comp_sol,
            I = trans.I_comp_sol,
            E = trans.E_comp_sol,
            CH = trans.CH_comp_sol,
            H_inc_comp_sol = trans.H_inc_comp_sol,
            H_prev_means = trans.H_prev_means,
            rt_init = trans.Rt_init
        )


    end


## model w/ wastewater and non-time varying hosp - not implemented
## model w/out wastewater and time-varying hosp - not implemented