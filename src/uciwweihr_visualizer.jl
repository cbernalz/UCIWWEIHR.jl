"""
    uciwweihr_visualizer(...)
Default visualizer for results of the UCIWWEIHR model, includes posterior/priors of generated quantities and posterior predictive samples for forecasting.  Forecasting plots will have the observed data alongside.

# Arguments
- `build_params::uciwweihr_model_params`: A uciwweihr_model_params object, if user desires priors next to plot (do not specify if you do not want prior plots).
- `data_hosp`: An array of hospital data used for model fitting, if user desires priors next to plot (do not specify if you do not want prior plots).
- `data_wastewater`: An array of wastewater data used for model fitting, if model does not use this do not specify this, if user desires priors next to plot (do not specify if you do not want prior plots).
- `obstimes`: An array of timepoints for observed hosp/wastewater, if user desires priors next to plot (do not specify if you do not want prior plots).
- `param_change_times`: An array of timepoints where the parameters change, if user desires priors next to plot (do not specify if you do not want prior plots).
- `seed`: Seed for the random number generator, if user desires priors next to plot (do not specify if you do not want prior plots).
- `forecast`: A boolean to indicate if forecasting is to be done, if user desires priors next to plot (do not specify if you do not want prior plots).
- `forecast_weeks`: Number of weeks to forecast, if user desires priors next to plot (do not specify if you do not want prior plots).
- `pp_sampeles`: Posterior predictive samples from the posterior/prior distribution, index 1 in uciwweihr_gq_pp output.
- `gq_samples`: Generated quantities samples from the posterior/prior distribution, index 2 in uciwweihr_gq_pp output.
- `obs_data_hosp`: An array of hospital data, data used for model fitting or extened timeseries for evaluation of forecast.
- `obs_data_wastewater`: An array of wastewater data, data used for model fitting or extened timeseries for evaluation of forecast.
- `actual_rt_vals`: An array of actual Rt values if user has access to them assumed to be on a daily scale.  This typically will come from some simulation.  Default is nothing.
- `actual_w_t`: An array of actual w_t values if user has access to them assumed to be on a daily scale.  This typically will come from some simulation.  Default is nothing.
- `actual_E_ode_sol`: An array of actual E ODE values if user has access to them assumed to be on a daily scale.  This typically will come from some simulation.  Default is nothing.
- `actual_I_ode_sol`: An array of actual I ODE values if user has access to them assumed to be on a daily scale.  This typically will come from some simulation.  Default is nothing.
- `actual_H_ode_sol`: An array of actual H ODE values if user has access to them assumed to be on a daily scale.  This typically will come from some simulation.  Default is nothing.
- `actual_non_time_varying_vals::uciwweihr_sim_params`: A uciwweihr_sim_params object of actual non-time varying parameter values if user has access to them. Default is nothing.
- `desired_params`: A list of lists of parameters to visualize. Each list will be visualized in a separate plot. Default is [["E_init", "I_init", "H_init"], ["gamma", "nu", "epsilon"], ["rho_gene", "tau", "df"], ["sigma_hosp"]].
- `time_varying_params`: A list of time varying parameters to visualize. Default is ["rt_vals", "w_t"].
- `var_to_pred`: A list of variables to predict. Default is ["data_wastewater", "data_hosp"].
- `quantiles`: A list of quantiles to calculate for ploting uncertainty. Default is [0.5, 0.8, 0.95].
- `bayes_dist_type`: A string to indicate if user is using Posterior or Prior distribution ("Posterior" / "Prior").
- `mcmcdaigs::Bool=true`: A boolean to indicate if user wants to visualize mcmc diagnosis plots and Effective Sample Size(ESS).
- `time_varying_plots::Bool=true`: A boolean to indicate if user wants to visualize time varying parameters.    
- `non_time_varying_plots::Bool=true`: A boolean to indicate if user wants to visualize non-time varying parameters.
- `ode_sol_plots::Bool=true`: A boolean to indicate if user wants to visualize ODE solutions.
- `pred_param_plots::Bool=true`: A boolean to indicate if user wants to visualize posterior (or prior) predictive parameter results.
- `save_plots::Bool=false`: A boolean to indicate if user wants to save the plots as pngs into a plots folder.
- `plot_name_to_save_mcmcdiag`: A string to indicate the name of the plot to save for MCMC diagnostics. Default is "mcmc_diagnosis_plots".
- `plot_name_to_save_time_varying`: A string to indicate the name of the plot to save for time varying parameters. Default is "mcmc_time_varying_parameter_plots".
- `plot_name_to_save_non_time_varying`: A string to indicate the name of the plot to save for non-time varying parameters. Default is "mcmc_nontime_varying_parameter_plots".
- `plot_name_to_save_ode_sol`: A string to indicate the name of the plot to save for ODE solutions. Default is "mcmc_ode_solution_plots".
- `plot_name_to_save_pred_param`: A string to indicate the name of the plot to save for posterior (or prior) predictive parameter results. Default is "mcmc_pred_parameter_plots".
"""
function uciwweihr_visualizer(
    data_hosp, 
    forecast_days,
    obstimes_hosp,
    param_change_times,
    seed,
    forecast,
    build_params::model_params_non_time_var_hosp_no_ww;
    pp_samples = nothing,
    gq_samples = nothing,
    samples = nothing,
    obs_data_hosp = nothing,
    obs_data_wastewater = nothing,
    actual_rt_vals = nothing,
    actual_w_t = nothing,
    actual_E_ode_sol = nothing,
    actual_I_ode_sol = nothing,
    actual_H_ode_sol = nothing,
    actual_non_time_varying_vals::uciwweihr_sim_params = uciwweihr_sim_params(ntuple(x->nothing, fieldcount(uciwweihr_sim_params))...),
    desired_params = [
        ["E_init", "I_init", "H_init"],
        ["gamma", "nu", "epsilon", "w"],
        ["rt_init", "sigma_Rt"],
        ["sigma_hosp"]
    ],
    time_varying_params = ["rt_vals"],
    var_to_pred = ["data_hosp"],
    quantiles = [0.5, 0.8, 0.95],
    bayes_dist_type = nothing,
    mcmcdaigs::Bool = true,
    time_varying_plots::Bool = true,
    non_time_varying_plots::Bool = true,
    ode_sol_plots::Bool = true,
    pred_param_plots::Bool = true,
    log_like_plots::Bool = true,
    save_plots::Bool = false,
    plot_name_to_save_mcmcdiag = "plots/mcmc_diagnosis_plots",
    plot_name_to_save_time_varying = "plots/mcmc_time_varying_parameter_plots",
    plot_name_to_save_non_time_varying = "plots/mcmc_nontime_varying_parameter_plots",
    plot_name_to_save_ode_sol = "plots/mcmc_ode_solution_plots",
    plot_name_to_save_pred_param = "plots/mcmc_pred_parameter_plots",
    plot_name_to_save_log_like = "plots/mcmc_log_likelihood_plots"
    )
    # Visualizer without wastewater data and constant hospitalization probability
    # Posterior/Prior Samples
    ## MCMC evaluation
    if mcmcdaigs
        mcmcdiags_vis(
            gq_samples, 
            desired_params, 
            actual_non_time_varying_vals;
            save_plots = save_plots,
            plot_name_to_save = plot_name_to_save_mcmcdiag
        )
    else
        println("MCMC Diagnostics Plots are not requested.")
    end

    if time_varying_plots
        time_varying_param_vis(
            build_params,
            data_hosp, 
            obstimes_hosp,
            param_change_times,
            seed,
            forecast,
            forecast_days;
            gq_samples = gq_samples,
            actual_rt_vals = actual_rt_vals,
            actual_w_t = actual_w_t,
            time_varying_params = time_varying_params,
            quantiles = quantiles,
            save_plots = save_plots,
            plot_name_to_save = plot_name_to_save_time_varying
        )
    else
        println("MCMC time varying parameter results are not requested.")
    end

    if non_time_varying_plots
        non_time_varying_param_vis(
            build_params,
            data_hosp, 
            obstimes_hosp,
            param_change_times,
            seed,
            forecast,
            forecast_days;
            gq_samples = gq_samples,
            desired_params = desired_params,
            actual_non_time_varying_vals = actual_non_time_varying_vals,
            save_plots = save_plots,
            plot_name_to_save = plot_name_to_save_non_time_varying
        )
    else
        println("MCMC non-time varying parameter results are not requested.")
    end

    if ode_sol_plots
        ode_solution_vis(
            gq_samples = gq_samples,
            actual_E_ode_sol = actual_E_ode_sol,
            actual_I_ode_sol = actual_I_ode_sol,
            actual_H_ode_sol = actual_H_ode_sol,
            quantiles = quantiles,
            save_plots = save_plots,
            plot_name_to_save = plot_name_to_save_ode_sol
        )
    else
        println("ODE Solution Plots are not requested.")
    end

    if pred_param_plots
        predictive_param_vis(
            pp_samples = pp_samples,
            data_wastewater = obs_data_wastewater,
            data_hosp = obs_data_hosp,
            obstimes_wastewater = obstimes_wastewater,
            obstimes_hosp = obstimes_hosp,
            forecast_days = forecast_days,
            vars_to_pred = var_to_pred, 
            quantiles = quantiles,
            bayes_dist_type = bayes_dist_type,
            save_plots = save_plots,
            plot_name_to_save = plot_name_to_save_pred_param
        )
    else
        println("MCMC posterior (or prior) predictive parameter results are not requested.")
    end

    if log_like_plots
        mcmcdiags_vis(
            samples, 
            [["lp"]];
            save_plots = save_plots,
            plot_name_to_save = plot_name_to_save_log_like
        )
    else
        println("MCMC Log Probs trace plot are not requested.")
    end
end


function uciwweihr_visualizer(
    data_hosp, 
    data_wastewater,
    forecast_days,
    obstimes_hosp,
    obstimes_wastewater,
    param_change_times,
    seed,
    forecast,
    build_params::model_params_time_var_hosp;
    pp_samples = nothing,
    gq_samples = nothing,
    samples = nothing,
    obs_data_hosp = nothing,
    obs_data_wastewater = nothing,
    actual_rt_vals = nothing,
    actual_w_t = nothing,
    actual_E_ode_sol = nothing,
    actual_I_ode_sol = nothing,
    actual_H_ode_sol = nothing,
    actual_non_time_varying_vals::uciwweihr_sim_params = uciwweihr_sim_params(ntuple(x->nothing, fieldcount(uciwweihr_sim_params))...),
    desired_params = [
        ["E_init", "I_init", "H_init"],
        ["gamma", "nu", "epsilon"],
        ["rt_init", "w_init"],
        ["sigma_w", "sigma_Rt"],
        ["sigma_ww", "sigma_hosp"], 
        ["rho_gene"]
    ],
    time_varying_params = ["rt_vals", "w_t"],
    var_to_pred = ["data_wastewater", "data_hosp"],
    quantiles = [0.5, 0.8, 0.95],
    bayes_dist_type = nothing,
    mcmcdaigs::Bool = true,
    time_varying_plots::Bool = true,
    non_time_varying_plots::Bool = true,
    ode_sol_plots::Bool = true,
    pred_param_plots::Bool = true,
    log_like_plots::Bool = true,
    save_plots::Bool = false,
    plot_name_to_save_mcmcdiag = "plots/mcmc_diagnosis_plots",
    plot_name_to_save_time_varying = "plots/mcmc_time_varying_parameter_plots",
    plot_name_to_save_non_time_varying = "plots/mcmc_nontime_varying_parameter_plots",
    plot_name_to_save_ode_sol = "plots/mcmc_ode_solution_plots",
    plot_name_to_save_pred_param = "plots/mcmc_pred_parameter_plots",
    plot_name_to_save_log_like = "plots/mcmc_log_likelihood_plots"
    )

    # Posterior/Prior Samples
    ## MCMC evaluation
    if mcmcdaigs
        mcmcdiags_vis(
            gq_samples, 
            desired_params, 
            actual_non_time_varying_vals;
            save_plots = save_plots,
            plot_name_to_save = plot_name_to_save_mcmcdiag
        )
    else
        println("MCMC Diagnostics Plots are not requested.")
    end

    if time_varying_plots
        time_varying_param_vis(
            build_params,
            data_hosp, 
            data_wastewater,
            obstimes_hosp,
            obstimes_wastewater,
            param_change_times,
            seed,
            forecast,
            forecast_days;
            gq_samples = gq_samples,
            actual_rt_vals = actual_rt_vals,
            actual_w_t = actual_w_t,
            time_varying_params = time_varying_params,
            quantiles = quantiles,
            save_plots = save_plots,
            plot_name_to_save = plot_name_to_save_time_varying
        )
    else
        println("MCMC time varying parameter results are not requested.")
    end

    if non_time_varying_plots
        non_time_varying_param_vis(
            build_params,
            data_hosp, 
            data_wastewater,
            obstimes_hosp,
            obstimes_wastewater,
            param_change_times,
            seed,
            forecast,
            forecast_days;
            gq_samples = gq_samples,
            desired_params = desired_params,
            actual_non_time_varying_vals = actual_non_time_varying_vals,
            save_plots = save_plots,
            plot_name_to_save = plot_name_to_save_non_time_varying
        )
    else
        println("MCMC non-time varying parameter results are not requested.")
    end

    if ode_sol_plots
        ode_solution_vis(
            gq_samples = gq_samples,
            actual_E_ode_sol = actual_E_ode_sol,
            actual_I_ode_sol = actual_I_ode_sol,
            actual_H_ode_sol = actual_H_ode_sol,
            quantiles = quantiles,
            save_plots = save_plots,
            plot_name_to_save = plot_name_to_save_ode_sol
        )
    else
        println("ODE Solution Plots are not requested.")
    end

    if pred_param_plots
        predictive_param_vis(
            pp_samples = pp_samples,
            data_wastewater = obs_data_wastewater,
            data_hosp = obs_data_hosp,
            obstimes_wastewater = obstimes_wastewater,
            obstimes_hosp = obstimes_hosp,
            forecast_days = forecast_days,
            vars_to_pred = var_to_pred, 
            quantiles = quantiles,
            bayes_dist_type = bayes_dist_type,
            save_plots = save_plots,
            plot_name_to_save = plot_name_to_save_pred_param
        )
    else
        println("MCMC posterior (or prior) predictive parameter results are not requested.")
    end

    if log_like_plots
        mcmcdiags_vis(
            samples, 
            [["lp"]];
            save_plots = save_plots,
            plot_name_to_save = plot_name_to_save_log_like
        )
    else
        println("MCMC Log Probs trace plot are not requested.")
    end
end


## model with wastewater and without time-varying hospitalization probability - not implemented
## model without wastewaeter and with time-varying hospitalization probability - not implemented


function uciwweihr_visualizer(
    forecast_days;
    pp_samples = nothing,
    gq_samples = nothing,
    samples = nothing,
    obs_data_hosp = nothing,
    obs_data_wastewater = nothing,
    obstimes_hosp = nothing,
    obstimes_wastewater = nothing,
    actual_rt_vals = nothing,
    actual_w_t = nothing,
    actual_E_ode_sol = nothing,
    actual_I_ode_sol = nothing,
    actual_H_ode_sol = nothing,
    actual_non_time_varying_vals::uciwweihr_sim_params = uciwweihr_sim_params(ntuple(x->nothing, fieldcount(uciwweihr_sim_params))...),
    desired_params = [
        ["E_init", "I_init", "H_init"],
        ["gamma", "nu", "epsilon"],
        ["rt_init", "w_init"],
        ["sigma_w", "sigma_Rt"],
        ["rho_gene"]
    ],
    time_varying_params = ["rt_vals", "w_t"],
    var_to_pred = ["data_wastewater", "data_hosp"],
    quantiles = [0.5, 0.8, 0.95],
    bayes_dist_type = nothing,
    mcmcdaigs::Bool = true,
    time_varying_plots::Bool = true,
    non_time_varying_plots::Bool = true,
    ode_sol_plots::Bool = true,
    pred_param_plots::Bool = true,
    save_plots::Bool = false,
    plot_name_to_save_mcmcdiag = "plots/mcmc_diagnosis_plots",
    plot_name_to_save_time_varying = "plots/mcmc_time_varying_parameter_plots",
    plot_name_to_save_non_time_varying = "plots/mcmc_nontime_varying_parameter_plots",
    plot_name_to_save_ode_sol = "plots/mcmc_ode_solution_plots",
    plot_name_to_save_pred_param = "plots/mcmc_pred_parameter_plots",
    plot_name_to_save_log_like = "plots/mcmc_log_likelihood_plots"
    )

    # Posterior/Prior Samples
    ## MCMC evaluation
    if mcmcdaigs
        mcmcdiags_vis(
            gq_samples, 
            desired_params, 
            actual_non_time_varying_vals;
            save_plots = save_plots,
            plot_name_to_save = plot_name_to_save_mcmcdiag
        )
    else
        println("MCMC Diagnostics Plots are not requested.")
    end

    if time_varying_plots
        time_varying_param_vis(
            gq_samples = gq_samples,
            actual_rt_vals = actual_rt_vals,
            actual_w_t = actual_w_t,
            time_varying_params = time_varying_params,
            quantiles = quantiles,
            save_plots = save_plots,
            plot_name_to_save = plot_name_to_save_time_varying
        )
    else
        println("MCMC time varying parameter results are not requested.")
    end

    if non_time_varying_plots
        non_time_varying_param_vis(
            gq_samples = gq_samples,
            desired_params = desired_params,
            actual_non_time_varying_vals = actual_non_time_varying_vals,
            save_plots = save_plots,
            plot_name_to_save = plot_name_to_save_non_time_varying
        )
    else
        println("MCMC non-time varying parameter results are not requested.")
    end

    if ode_sol_plots
        ode_solution_vis(
            gq_samples = gq_samples,
            actual_E_ode_sol = actual_E_ode_sol,
            actual_I_ode_sol = actual_I_ode_sol,
            actual_H_ode_sol = actual_H_ode_sol,
            quantiles = quantiles,
            save_plots = save_plots,
            plot_name_to_save = plot_name_to_save_ode_sol
        )
    else
        println("ODE Solution Plots are not requested.")
    end

    if pred_param_plots
        predictive_param_vis(
            pp_samples = pp_samples,
            data_wastewater = obs_data_wastewater,
            data_hosp = obs_data_hosp,
            obstimes_wastewater = obstimes_wastewater,
            obstimes_hosp = obstimes_hosp,
            forecast_days = forecast_days,
            vars_to_pred = var_to_pred, 
            quantiles = quantiles,
            bayes_dist_type = bayes_dist_type,
            save_plots = save_plots,
            plot_name_to_save = plot_name_to_save_pred_param
        )
    else
        println("MCMC posterior (or prior) predictive parameter results are not requested.")
    end

    if log_like_plots
        mcmcdiags_vis(
            samples, 
            [["lp"]];
            save_plots = save_plots,
            plot_name_to_save = plot_name_to_save_log_like
        )
    else
        println("MCMC Log Probs trace plot are not requested.")
    end
end
