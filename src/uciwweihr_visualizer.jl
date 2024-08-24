"""
    uciwweihr_visualizer(...)
Default visualizer for results of the UCIWWEIHR model, includes posterior/priors of generated quantities and posterior predictive samples for forecasting.  Forecasting plots will have the observed data alongside.

# Arguments
- `pp_sampeles`: Posterior predictive samples from the posterior/prior distribution, index 1 in uciwweihr_gq_pp output.
- `gq_samples`: Generated quantities samples from the posterior/prior distribution, index 2 in uciwweihr_gq_pp output.
- `data_hosp`: An array of hospital data.
- `data_wastewater`: An array of pathogen genome concentration in localized wastewater data.
- `actual_rt_vals`: An array of actual Rt values if user has access to them assumed to be on a daily scale.  This typically will come from some simulation.  Default is nothing.
- `actual_w_t`: An array of actual w_t values if user has access to them assumed to be on a daily scale.  This typically will come from some simulation.  Default is nothing.
- `actual_non_time_varying_vals::uciwweihr_sim_params`: A uciwweihr_sim_params object of actual non-time varying parameter values if user has access to them. Default is nothing.
- `forecast_weeks`: Number of weeks to forecasted.  Default is 0.
- `desired_params`: A list of lists of parameters to visualize. Each list will be visualized in a separate plot. Default is [["E_init", "I_init", "H_init"], ["gamma", "nu", "epsilon"], ["rho_gene", "tau", "df"], ["sigma_hosp"]].
- `time_varying_params`: A list of time varying parameters to visualize. Default is ["rt_vals", "w_t"].
- `var_to_pred`: A list of variables to predict. Default is ["data_wastewater", "data_hosp"].
- `quantiles`: A list of quantiles to calculate for ploting uncertainty. Default is [0.5, 0.8, 0.95].
- `bayes_dist_type`: A string to indicate if user is using Posterior or Prior distribution ("Posterior" / "Prior").
- `mcmcdaigs::Bool=true`: A boolean to indicate if user wants to visualize mcmc diagnosis plots and Effective Sample Size(ESS).
- `time_varying_plots::Bool=true`: A boolean to indicate if user wants to visualize time varying parameters.    
- `non_time_varying_plots::Bool=true`: A boolean to indicate if user wants to visualize non-time varying parameters.
- `pred_param_plots::Bool=true`: A boolean to indicate if user wants to visualize posterior (or prior) predictive parameter results.
- `save_plots::Bool=false`: A boolean to indicate if user wants to save the plots as pngs into a plots folder.
"""
function uciwweihr_visualizer(;
    pp_samples = nothing,
    gq_samples = nothing,
    data_hosp = nothing,
    data_wastewater = nothing,
    actual_rt_vals = nothing,
    actual_w_t = nothing,
    actual_non_time_varying_vals::uciwweihr_sim_params = nothing,
    forecast_weeks = 0,
    desired_params = [
        ["E_init", "I_init", "H_init"],
        ["gamma", "nu", "epsilon"],
        ["rt_init", "w_init"],
        ["rho_gene", "tau", "df"],
        ["sigma_hosp"]
    ],
    time_varying_params = ["rt_vals", "w_t"],
    var_to_pred = ["data_wastewater", "data_hosp"],
    quantiles = [0.5, 0.8, 0.95],
    bayes_dist_type = nothing,
    mcmcdaigs::Bool = true,
    time_varying_plots::Bool = true,
    non_time_varying_plots::Bool = true,
    pred_param_plots::Bool = true,
    save_plots::Bool = false
    )

    # Posterior/Prior Samples
    ## MCMC evaluation
    if mcmcdaigs
        mcmcdiags_vis(
            gq_samples = gq_samples, 
            desired_params = desired_params, 
            actual_non_time_varying_vals = actual_non_time_varying_vals,
            save_plots = save_plots
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
            save_plots = save_plots
        )
    else
        println("MCMC time varying parameter results are not requested.")
    end

    if non_time_varying_plots
        non_time_varying_param_vis(
            gq_samples = gq_samples,
            desired_params = desired_params,
            bayes_dist_type = bayes_dist_type,
            actual_non_time_varying_vals = actual_non_time_varying_vals,
            save_plots = save_plots
        )
    else
        println("MCMC non-time varying parameter results are not requested.")
    end

    if pred_param_plots
        predictive_param_vis(
            pp_samples = pp_samples,
            data_wastewater = data_wastewater,
            data_hosp = data_hosp,
            forecast_weeks = forecast_weeks,
            vars_to_pred = var_to_pred, 
            quantiles = quantiles,
            bayes_dist_type = bayes_dist_type,
            save_plots = save_plots
        )
    else
        println("MCMC posterior (or prior) predictive parameter results are not requested.")
    end



end
