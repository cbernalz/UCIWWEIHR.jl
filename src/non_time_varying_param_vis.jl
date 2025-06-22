"""
    non_time_varying_param_vis(...)

Used in the `uciwweihr_visualizer` to create visuals for non-time varying parameters.

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
- `desired_params`: A list of lists of parameters to visualize. Each list will be visualized in a separate plot. Default is any parameter not in this list : ["alpha_t", "w_t", "rt_vals", "log_genes_mean", "H"]
- `actual_non_time_varying_vals::uciwweihr_sim_params`: A uciwweihr_sim_params object of actual non-time varying parameter values if user has access to them. Default is nothing.
- `save_plots::Bool=false`: A boolean to indicate if user wants to save the plots as pngs into a plots folder.
- `plot_name_to_save`: A string to indicate the name of the plot to save. Default is "mcmc_nontime_varying_parameter_plots".
"""
function non_time_varying_param_vis(
    build_params::model_params_non_time_var_hosp_no_ww,
    data_hosp, 
    obstimes_hosp,
    param_change_times,
    seed,
    forecast,
    forecast_days;
    gq_samples=nothing,
    desired_params=nothing,
    actual_non_time_varying_vals::uciwweihr_sim_params = nothing,
    save_plots::Bool=false,
    plot_name_to_save = "mcmc_nontime_varying_parameter_plots"
    )
    println("Generating non-time varying parameter plots (w/out wastewater, w const hospitalization probability)...")
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

    non_time_varying_plots = []    
    for param_group in desired_params
        for curr_param in param_group
            chains = unique(gq_samples.chain)
            for chain in chains
                curr_param_chain_df = filter(row -> row.chain == chain, gq_samples)
                prior_param_samples = filter(row -> row.chain == chain, gq_samples_priors)
                plt = histogram(prior_param_samples[:, curr_param], 
                           label = "Chain $chain (Prior)", 
                           bins = 50,
                            normalize = :probability,
                           color = :blue,
                           alpha = 0.1)
                histogram!(plt, curr_param_chain_df[:, curr_param], 
                                label = "Chain $chain (Posterior)", 
                                title = "$curr_param",
                                bins = 50,
                                normalize = :probability,
                                xlabel = "Value for $curr_param",
                                ylabel = "Probability",
                                alpha = 0.7,
                                color = :blue,
                                legend = :topright)
                if !isnothing(actual_non_time_varying_vals.time_points)
                    actual_param_value = round(getfield(actual_non_time_varying_vals, Symbol(curr_param)), digits=3)
                    vline!(plt, 
                            [actual_param_value],
                            label = "Actual Value: $actual_param_value",
                            color = :red, 
                            linewidth = 3,
                            linestyle=:dash)
                end
                
                push!(non_time_varying_plots, plt) 
            end
        end
    end
    
    if !isempty(non_time_varying_plots)
        num_plots = length(non_time_varying_plots)
        num_chains = length(unique(gq_samples.chain))
        num_params = length(desired_params)
        layout_rows = num_chains * length(desired_params[1])
        layout_cols = num_params
    
        final_plot = plot(non_time_varying_plots..., 
                            layout = (layout_rows, layout_cols), 
                            size = (1500, 1500))
        display(final_plot)
        if save_plots
            save_plots_to_docs(final_plot, plot_name_to_save)
        end
    else
        println("NO NON-TIME VARYING PARAMETER PLOTS TO DISPLAY!!!")
    end

end


function non_time_varying_param_vis(
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
    desired_params=nothing,
    actual_non_time_varying_vals::uciwweihr_sim_params = nothing,
    save_plots::Bool=false,
    plot_name_to_save = "mcmc_nontime_varying_parameter_plots"
    )
    println("Generating non-time varying parameter plots (with wastewater and with time-varying hospitalization probability)...")
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

    non_time_varying_plots = []    
    for param_group in desired_params
        for curr_param in param_group
            chains = unique(gq_samples.chain)
            for chain in chains
                curr_param_chain_df = filter(row -> row.chain == chain, gq_samples)
                prior_param_samples = filter(row -> row.chain == chain, gq_samples_priors)
                plt = histogram(prior_param_samples[:, curr_param], 
                           label = "Chain $chain (Prior)", 
                           bins = 50,
                            normalize = :probability,
                           color = :blue,
                           alpha = 0.1)
                histogram!(plt, curr_param_chain_df[:, curr_param], 
                                label = "Chain $chain (Posterior)", 
                                title = "$curr_param",
                                bins = 50,
                                normalize = :probability,
                                xlabel = "Value for $curr_param",
                                ylabel = "Probability",
                                alpha = 0.7,
                                color = :blue,
                                legend = :topright)
                if !isnothing(actual_non_time_varying_vals.time_points)
                    actual_param_value = round(getfield(actual_non_time_varying_vals, Symbol(curr_param)), digits=3)
                    vline!(plt, 
                            [actual_param_value],
                            label = "Actual Value: $actual_param_value",
                            color = :red, 
                            linewidth = 3,
                            linestyle=:dash)
                end
                
                push!(non_time_varying_plots, plt) 
            end
        end
    end
    
    if !isempty(non_time_varying_plots)
        num_plots = length(non_time_varying_plots)
        num_chains = length(unique(gq_samples.chain))
        num_params = length(desired_params)
        layout_rows = num_chains * length(desired_params[1])
        layout_cols = num_params
    
        final_plot = plot(non_time_varying_plots..., 
                            layout = (layout_rows, layout_cols), 
                            size = (1500, 1500))
        display(final_plot)
        if save_plots
            save_plots_to_docs(final_plot, plot_name_to_save)
        end
    else
        println("NO NON-TIME VARYING PARAMETER PLOTS TO DISPLAY!!!")
    end

end

## model with wastewater and without time-varying hospitalization probability - not implemented
## model without wastewater and with time-varying hospitalization probability - not implemented


function non_time_varying_param_vis(;
    gq_samples=nothing,
    desired_params=nothing,
    actual_non_time_varying_vals::uciwweihr_sim_params = nothing,
    save_plots::Bool=false,
    plot_name_to_save = "mcmc_nontime_varying_parameter_plots"
    )
    # Plotting non time varying parameters
    non_time_varying_plots = []
    for param_group in desired_params
        for curr_param in param_group
            chains = unique(gq_samples.chain)
            for chain in chains
                curr_param_chain_df = filter(row -> row.chain == chain, gq_samples)
                plt = histogram(curr_param_chain_df[:, curr_param], 
                                label = "Chain $chain", 
                                title = "$curr_param",
                                bins = 50,
                                normalize = :probability,
                                xlabel = "Probability",
                                ylabel = "Value for $curr_param",
                                xguidefont = font(8),
                                yguidefont = font(8), 
                                titlefont = font(10), 
                                legendfont = font(8),
                                legend = :topright)
                if !isnothing(actual_non_time_varying_vals.time_points)
                    actual_param_value = round(getfield(actual_non_time_varying_vals, Symbol(curr_param)), digits=3)
                    vline!(plt, 
                            [actual_param_value],
                            label = "Actual Value: $actual_param_value",
                            color = :red, 
                            linewidth = 3,
                            linestyle=:dash)
                end
                push!(non_time_varying_plots, plt) 
            end
        end
    end
    if !isempty(non_time_varying_plots)
        num_plots = length(non_time_varying_plots)
        num_chains = length(unique(gq_samples.chain))
        num_params = length(desired_params)
        layout_rows = num_chains * length(desired_params[1])
        layout_cols = num_params
    
        final_plot = plot(non_time_varying_plots..., 
                            layout = (layout_rows, layout_cols), 
                            size = (1500, 1500))
        display(final_plot)
        if save_plots
            save_plots_to_docs(final_plot, plot_name_to_save)
        end
    else
        println("NO NON-TIME VARYING PARAMETER PLOTS TO DISPLAY!!!")
    end
    
    
end

