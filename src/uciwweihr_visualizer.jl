"""
    uciwweihr_visualizer(...)
Default visualizer for results of the UCIWWEIHR model, includes posterior/priors of generated quantities and posterior predictive samples for forecasting.  Forecasting plots will have the observed data alongside.

# Arguments
- `pp_sampeles`: Posterior predictive samples from the posterior/prior distribution, index 1 in uciwweihr_gq_pp output.
- `gq_samples`: Generated quantities samples from the posterior/prior distribution, index 2 in uciwweihr_gq_pp output.
- `data_hosp`: An array of hospital data.
- `data_wastewater`: An array of pathogen genome concentration in localized wastewater data.
- `obstimes`: An array of timepoints for observed hosp/wastewater.
- `param_change_times`: An array of timepoints where the parameters change.
- `actual_rt_vals`: An array of actual Rt values if user has access to them assumed to be on a daily scale.  This typically will come from some simulation.  Default is nothing.
- `actual_w_t`: An array of actual w_t values if user has access to them assumed to be on a daily scale.  This typically will come from some simulation.  Default is nothing.
- `desired_params`: A list of lists of parameters to visualize. Each list will be visualized in a separate plot. Default is [["E_init", "I_init", "H_init"], ["gamma", "nu", "epsilon"], ["rho_gene", "tau", "df"], ["sigma_hosp"]].
- `time_varying_params`: A list of time varying parameters to visualize. Default is ["rt_vals", "w_t"].
- `quantiles`: A list of quantiles to calculate for ploting uncertainty. Default is [0.5, 0.8, 0.95].
- `mcmcdaigs::Bool=true`: A boolean to indicate if user wants to visualize mcmc diagnosis plots and Effective Sample Size(ESS).
- `time_varying_plots::Bool=true`: A boolean to indicate if user wants to visualize time varying parameters.    
- `save_plots::Bool=false`: A boolean to indicate if user wants to save the plots as pngs into a plots folder.
"""
function uciwweihr_visualizer(;
    pp_samples=nothing,
    gq_samples=nothing,
    data_hosp=nothing,
    data_wastewater=nothing,
    obstimes=nothing,
    param_change_times=nothing,
    actual_rt_vals=nothing,
    actual_w_t=nothing,
    desired_params=[
        ["E_init", "I_init", "H_init"],
        ["gamma", "nu", "epsilon"],
        ["rt_init", "w_init"],
        ["rho_gene", "tau", "df"],
        ["sigma_hosp"]
    ],
    time_varying_params = ["rt_vals", "w_t"],
    quantiles = [0.5, 0.8, 0.95],
    mcmcdaigs::Bool=true,
    time_varying_plots::Bool=true,
    save_plots::Bool=false
    )

    # Posterior/Prior Samples
    ## MCMC evaluation
    if mcmcdaigs
        mcmcdiags_vis(
            gq_samples = gq_samples, 
            desired_params = desired_params, 
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



end
