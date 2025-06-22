"""
    time_varying_param_vis(...)
Used in the `uciwweihr_visualizer` to create visuals for time varying parameters.

# Arguments
- `build_params::uciwweihr_model_params`: A struct of model parameters used to build `gq_samples`, used only if user desired priors next to posteriors.
- `data_hosp`: Hospitalization data, used only if user desired priors next to posteriors.
- `data_wastewater`: Wastewater data, if model does not use this do not specify this, if user desires priors next to plot (do not specify if you do not want prior plots).
- `obstimes_hosp`: An array of time points for hospital data, used only if user desired priors next to posteriors.
- `obstimes_wastewater`: An array of time points for wastewater data, used only if user desired priors next to posteriors.
- `param_change_times`: An array of time points where the parameters change, used only if user desired priors next to posteriors.
- `seed`: An integer to set the seed for reproducibility, used only if user desired priors next to posteriors.
- `forecast`: A boolean to indicate if user wants to forecast, used only if user desired priors next to posteriors.
- `forecast_weeks`: An integer to indicate the number of weeks to forecast, used only if user desired priors next to posteriors.
- `gq_samples`: Generated quantities samples from the posterior/prior distribution, index 2 in uciwweihr_gq_pp output.
- `actual_rt_vals`: An array of actual Rt values if user has access to them assumed to be on a daily scale.  This typically will come from some simulation.  Default is nothing.
- `actual_w_t`: An array of actual w_t values if user has access to them assumed to be on a daily scale.  This typically will come from some simulation.  Default is nothing.
- `time_varying_params`: A list of time varying parameters to visualize. Default is ["rt_vals", "w_t"].
- `quantiles`: A list of quantiles to calculate for ploting uncertainty. Default is [0.5, 0.8, 0.95].
- `save_plots::Bool=false`: A boolean to indicate if user wants to save the plots as pngs into a plots folder.
- `plot_name_to_save`: A string to indicate the name of the plot to save. Default is "mcmc_time_varying_parameter_plots".
"""

function time_varying_param_vis(
    build_params::model_params_non_time_var_hosp_no_ww,
    data_hosp, 
    obstimes_hosp,
    param_change_times,
    seed,
    forecast,
    forecast_days;
    gq_samples=nothing,
    actual_rt_vals=nothing,
    actual_w_t=nothing,
    time_varying_params = ["rt_vals"],
    quantiles = [0.5, 0.8, 0.95],
    save_plots::Bool=false,
    plot_name_to_save = "mcmc_time_varying_parameter_plots"
    )
    println("Generating time varying parameter plots (with priors, w/out wastewater, w const hospitalization probability)...")
    chains = unique(gq_samples.chain)
    samples = fit(
        data_hosp,
        obstimes_hosp,
        param_change_times,
        build_params;
        priors_only=true,
        n_chains=length(chains),
        seed=seed
        )
    prior_model_output = generate_pq_pp(
        samples,
        data_hosp,
        obstimes_hosp,
        param_change_times,
        build_params;
        seed=seed,
        forecast=forecast, 
        forecast_days=forecast_days
    )
    gq_samples_priors = prior_model_output[2]

    # Plotting time-varying parameters for both posterior and prior
    time_varying_plots = []
    column_names = names(gq_samples)
    column_names_priors = names(gq_samples_priors)

    for var_prefix in time_varying_params
        time_varying_param = filter(name -> startswith_any(name, [var_prefix]), column_names)
        time_varying_param_priors = filter(name -> startswith_any(name, [var_prefix]), column_names_priors)

        time_varying_subset_df = gq_samples[:, [time_varying_param..., "iteration", "chain"]]
        time_varying_subset_df_priors = gq_samples_priors[:, [time_varying_param_priors..., "iteration", "chain"]]

        chains = unique(time_varying_subset_df.chain)

        for chain in chains
            medians, lower_bounds, upper_bounds = calculate_quantiles(time_varying_subset_df, chain, var_prefix, quantiles)
            medians_priors, lower_bounds_priors, upper_bounds_priors = calculate_quantiles(time_varying_subset_df_priors, chain, var_prefix, quantiles)

            ribbon_colors = generate_colors(length(quantiles))
            daily_medians = repeat(medians, inner=7)
            daily_lower_bounds = repeat(lower_bounds, inner=7)
            daily_upper_bounds = repeat(upper_bounds, inner=7)

            daily_medians_priors = repeat(medians_priors, inner=7)
            daily_lower_bounds_priors = repeat(lower_bounds_priors, inner=7)
            daily_upper_bounds_priors = repeat(upper_bounds_priors, inner=7)

            daily_x = 1:length(daily_medians)
            plt = plot(title = "Quantiles for Chain $chain for $var_prefix (Posterior)",
                       xlabel = "Time Points (daily scale)",
                       ylabel = "Value for $var_prefix")

            for (i, q) in enumerate(quantiles)
                daily_upper_bounds_temp = map(x -> x[i], daily_upper_bounds)
                daily_lower_bounds_temp = map(x -> x[i], daily_lower_bounds)
                plot!(plt, daily_x, daily_medians, ribbon = (daily_upper_bounds_temp .- daily_medians, daily_medians .- daily_lower_bounds_temp),
                      fillalpha = 0.2,
                      label = "$(@sprintf("%.0f", q*100))% Quantile (Posterior)",
                      color = ribbon_colors[i],
                      fillcolor = ribbon_colors[i],
                      legend = :topright)
            end
            plot!(plt, daily_x, daily_medians, label = "Median (Posterior)", color = :black, lw = 2)

            if !isnothing(actual_rt_vals) && var_prefix == "rt_vals"
                scatter!(plt, 1:length(actual_rt_vals), actual_rt_vals, label = "Actual Rt Values", color = :red, lw = 2, marker = :circle)
            end
            if !isnothing(actual_w_t) && var_prefix == "w_t"
                scatter!(plt, 1:length(actual_w_t), actual_w_t, label = "Actual w_t Values", color = :red, lw = 2, marker = :circle)
            end

            plt_priors = plot(title = "Quantiles for Chain $chain for $var_prefix (Prior)",
                              xlabel = "Time Points (daily scale)",
                              ylabel = "Value for $var_prefix")

            for (i, q) in enumerate(quantiles)
                daily_upper_bounds_priors_temp = map(x -> x[i], daily_upper_bounds_priors)
                daily_lower_bounds_priors_temp = map(x -> x[i], daily_lower_bounds_priors)
                plot!(plt_priors, daily_x, daily_medians_priors, ribbon = (daily_upper_bounds_priors_temp .- daily_medians_priors, daily_medians_priors .- daily_lower_bounds_priors_temp),
                      fillalpha = 0.2,
                      label = "$(@sprintf("%.0f", q*100))% Quantile (Prior)",
                      color = ribbon_colors[i],
                      fillcolor = ribbon_colors[i],
                      legend = :topright)
            end
            plot!(plt_priors, daily_x, daily_medians_priors, label = "Median (Prior)", color = :blue, lw = 2)
            
            if !isnothing(actual_rt_vals) && var_prefix == "rt_vals"
                scatter!(plt_priors, 1:length(actual_rt_vals), actual_rt_vals, label = "Actual Rt Values", color = :red, lw = 2, marker = :circle)
            end
            if !isnothing(actual_w_t) && var_prefix == "w_t"
                scatter!(plt_priors, 1:length(actual_w_t), actual_w_t, label = "Actual w_t Values", color = :red, lw = 2, marker = :circle)
            end    


            push!(time_varying_plots, plt, plt_priors)
        end
    end

    if !isempty(time_varying_plots)
        chains = unique(gq_samples.chain)
        plt = plot(time_varying_plots..., layout = (length(time_varying_params), length(chains) * 2), size = (1000, 1000))
        display(plt)
        if save_plots
            save_plots_to_docs(plt, plot_name_to_save)
        end
    else
        println("NO TIME VARYING PARAMETER PLOTS TO DISPLAY!!!")
    end

    
end


function time_varying_param_vis(
    build_params::model_params_time_var_hosp,
    data_hosp, 
    data_wastewater,
    obstimes_hosp,
    obstimes_wastewater,
    param_change_times,
    seed,
    forecast,
    forecast_days;
    gq_samples=nothing,
    actual_rt_vals=nothing,
    actual_w_t=nothing,
    time_varying_params = ["rt_vals", "w_t"],
    quantiles = [0.5, 0.8, 0.95],
    save_plots::Bool=false,
    plot_name_to_save = "mcmc_time_varying_parameter_plots"
    )
    println("Generating time varying parameter plots (with wastewater and time-varying hospitalization probability)...")
    chains = unique(gq_samples.chain)
    samples = fit(
        data_hosp,
        data_wastewater,
        obstimes_hosp,
        obstimes_wastewater,
        param_change_times,
        build_params;
        priors_only=true,
        n_chains=length(chains),
        seed=seed
        )
    prior_model_output = generate_pq_pp(
        samples,
        data_hosp,
        data_wastewater,
        obstimes_hosp,
        obstimes_wastewater,
        param_change_times,
        build_params;
        seed=seed,
        forecast=forecast, forecast_days=forecast_days
    )
    gq_samples_priors = prior_model_output[2]

    # Plotting time-varying parameters for both posterior and prior
    time_varying_plots = []
    column_names = names(gq_samples)
    column_names_priors = names(gq_samples_priors)

    for var_prefix in time_varying_params
        time_varying_param = filter(name -> startswith_any(name, [var_prefix]), column_names)
        time_varying_param_priors = filter(name -> startswith_any(name, [var_prefix]), column_names_priors)

        time_varying_subset_df = gq_samples[:, [time_varying_param..., "iteration", "chain"]]
        time_varying_subset_df_priors = gq_samples_priors[:, [time_varying_param_priors..., "iteration", "chain"]]

        chains = unique(time_varying_subset_df.chain)

        for chain in chains
            medians, lower_bounds, upper_bounds = calculate_quantiles(time_varying_subset_df, chain, var_prefix, quantiles)
            medians_priors, lower_bounds_priors, upper_bounds_priors = calculate_quantiles(time_varying_subset_df_priors, chain, var_prefix, quantiles)

            ribbon_colors = generate_colors(length(quantiles))
            daily_medians = repeat(medians, inner=7)
            daily_lower_bounds = repeat(lower_bounds, inner=7)
            daily_upper_bounds = repeat(upper_bounds, inner=7)

            daily_medians_priors = repeat(medians_priors, inner=7)
            daily_lower_bounds_priors = repeat(lower_bounds_priors, inner=7)
            daily_upper_bounds_priors = repeat(upper_bounds_priors, inner=7)

            daily_x = 1:length(daily_medians)
            plt = plot(title = "Quantiles for Chain $chain for $var_prefix (Posterior)",
                       xlabel = "Time Points (daily scale)",
                       ylabel = "Value for $var_prefix")

            for (i, q) in enumerate(quantiles)
                daily_upper_bounds_temp = map(x -> x[i], daily_upper_bounds)
                daily_lower_bounds_temp = map(x -> x[i], daily_lower_bounds)
                plot!(plt, daily_x, daily_medians, ribbon = (daily_upper_bounds_temp .- daily_medians, daily_medians .- daily_lower_bounds_temp),
                      fillalpha = 0.2,
                      label = "$(@sprintf("%.0f", q*100))% Quantile (Posterior)",
                      color = ribbon_colors[i],
                      fillcolor = ribbon_colors[i],
                      legend = :topright)
            end
            plot!(plt, daily_x, daily_medians, label = "Median (Posterior)", color = :black, lw = 2)

            if !isnothing(actual_rt_vals) && var_prefix == "rt_vals"
                scatter!(plt, 1:length(actual_rt_vals), actual_rt_vals, label = "Actual Rt Values", color = :red, lw = 2, marker = :circle)
            end
            if !isnothing(actual_w_t) && var_prefix == "w_t"
                scatter!(plt, 1:length(actual_w_t), actual_w_t, label = "Actual w_t Values", color = :red, lw = 2, marker = :circle)
            end

            plt_priors = plot(title = "Quantiles for Chain $chain for $var_prefix (Prior)",
                              xlabel = "Time Points (daily scale)",
                              ylabel = "Value for $var_prefix")

            for (i, q) in enumerate(quantiles)
                daily_upper_bounds_priors_temp = map(x -> x[i], daily_upper_bounds_priors)
                daily_lower_bounds_priors_temp = map(x -> x[i], daily_lower_bounds_priors)
                plot!(plt_priors, daily_x, daily_medians_priors, ribbon = (daily_upper_bounds_priors_temp .- daily_medians_priors, daily_medians_priors .- daily_lower_bounds_priors_temp),
                      fillalpha = 0.2,
                      label = "$(@sprintf("%.0f", q*100))% Quantile (Prior)",
                      color = ribbon_colors[i],
                      fillcolor = ribbon_colors[i],
                      legend = :topright)
            end
            plot!(plt_priors, daily_x, daily_medians_priors, label = "Median (Prior)", color = :blue, lw = 2)
            
            if !isnothing(actual_rt_vals) && var_prefix == "rt_vals"
                scatter!(plt_priors, 1:length(actual_rt_vals), actual_rt_vals, label = "Actual Rt Values", color = :red, lw = 2, marker = :circle)
            end
            if !isnothing(actual_w_t) && var_prefix == "w_t"
                scatter!(plt_priors, 1:length(actual_w_t), actual_w_t, label = "Actual w_t Values", color = :red, lw = 2, marker = :circle)
            end    


            push!(time_varying_plots, plt, plt_priors)
        end
    end

    if !isempty(time_varying_plots)
        chains = unique(gq_samples.chain)
        plt = plot(time_varying_plots..., layout = (length(time_varying_params), length(chains) * 2), size = (1000, 1000))
        display(plt)
        if save_plots
            save_plots_to_docs(plt, plot_name_to_save)
        end
    else
        println("NO TIME VARYING PARAMETER PLOTS TO DISPLAY!!!")
    end

    
end

## model without wastewater and with time-varying hospitalization probability - not implemented
## model with wastewater and without time-varying hospitalization probability - not implemented

function time_varying_param_vis(;
    gq_samples=nothing,
    actual_rt_vals=nothing,
    actual_w_t=nothing,
    time_varying_params = ["rt_vals", "w_t"],
    quantiles = [0.5, 0.8, 0.95],
    save_plots::Bool=false,
    plot_name_to_save = "mcmc_time_varying_parameter_plots"
    )
    println("Generating time varying parameter plots (w/out priors)...")
    # Plotting time varying parameters
    time_varying_plots = []
    column_names = names(gq_samples)
    for var_prefix in time_varying_params
        time_varying_param = filter(name -> startswith_any(name, [var_prefix]), column_names)
        time_varying_subset_df = gq_samples[:, [time_varying_param..., "iteration", "chain"]]
        chains = unique(time_varying_subset_df.chain)
        for chain in chains
            medians, lower_bounds, upper_bounds = calculate_quantiles(time_varying_subset_df, chain, var_prefix, quantiles)
            ribbon_colors = generate_colors(length(quantiles))
            daily_medians = repeat(medians, inner=7)
            daily_lower_bounds = repeat(lower_bounds, inner=7)
            daily_upper_bounds = repeat(upper_bounds, inner=7)
            daily_x = 1:length(daily_medians)
            plt = plot(title = "Quantiles for Chain $chain for $var_prefix",
                        xlabel = "Time Points (daily scale)",
                        ylabel = "Value for $var_prefix")
            for (i, q) in enumerate(quantiles)
                daily_upper_bounds_temp = map(x -> x[i], daily_upper_bounds)
                daily_lower_bounds_temp = map(x -> x[i], daily_lower_bounds)
                plot!(plt, daily_x, daily_medians, ribbon = (daily_upper_bounds_temp .- daily_medians, daily_medians .- daily_lower_bounds_temp),
                        fillalpha = 0.2,
                        label = "$(@sprintf("%.0f", q*100))% Quantile",
                        color = ribbon_colors[i],
                       fillcolor = ribbon_colors[i])
            end
            plot!(plt, daily_x, daily_medians, label = "Median", color = :black, lw = 2)
            if !isnothing(actual_rt_vals) && var_prefix == "rt_vals"
                scatter!(plt, 1:length(actual_rt_vals), actual_rt_vals, label = "Actual Rt Values", color = :red, lw = 2, marker = :circle)
            end
            if !isnothing(actual_w_t) && var_prefix == "w_t"
                scatter!(plt, 1:length(actual_w_t), actual_w_t, label = "Actual w_t Values", color = :red, lw = 2, marker = :circle)
            end
            push!(time_varying_plots, plt) 
        end
    end

    if !isempty(time_varying_plots)
        chains = unique(gq_samples.chain)
        plt = plot(time_varying_plots..., layout = (length(time_varying_params), length(chains)), size = (1000, 1000))
        display(plt)
        if save_plots
            save_plots_to_docs(plt, plot_name_to_save)
        end
    else
        println("NO TIME VARYING PARAMETER PLOTS TO DISPLAY!!!")
    end
    
end
