"""
    non_time_varying_param_vis(...)

Used in the `uciwweihr_visualizer` to create visuals for non-time varying parameters.

# Arguments
- `gq_samples`: Generated quantities samples from the posterior/prior distribution, index 2 in uciwweihr_gq_pp output.
- `desired_params`: A list of lists of parameters to visualize. Each list will be visualized in a separate plot. Default is any parameter not in this list : ["alpha_t", "w_t", "rt_vals", "log_genes_mean", "H"]
- `bayes_dist_type`: A string to indicate if user is using Posterior or Prior distribution. Default is "Posterior".
- `actual_non_time_varying_vals::uciwweihr_sim_params`: A uciwweihr_sim_params object of actual non-time varying parameter values if user has access to them. Default is nothing.
- `save_plots::Bool=false`: A boolean to indicate if user wants to save the plots as pngs into a plots folder.
- `plot_name_to_save`: A string to indicate the name of the plot to save. Default is "mcmc_nontime_varying_parameter_plots".
"""
function non_time_varying_param_vis(
    build_params::uciwweihr_model_params,
    data_hosp, 
    obstimes,
    param_change_times,
    seed,
    forecast,
    forecast_weeks;
    gq_samples=nothing,
    desired_params=nothing,
    bayes_dist_type="Posterior",
    actual_non_time_varying_vals::uciwweihr_sim_params = nothing,
    save_plots::Bool=false,
    plot_name_to_save = "mcmc_nontime_varying_parameter_plots"
    )
    println("Generating non-time varying parameter plots (with priors and w/out wastewater)...")
    samples = uciwweihr_fit(
        data_hosp;
        obstimes = obstimes,
        param_change_times = param_change_times,
        priors_only = true,
        seed = seed,
        params = build_params
        )
    prior_model_output = uciwweihr_gq_pp(
        samples,
        data_hosp;
        obstimes,
        param_change_times,
        seed = seed,
        params = build_params,
        forecast = forecast, 
        forecast_weeks = forecast_weeks
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
                                title = "$bayes_dist_type $curr_param",
                                bins = 50,
                                normalize = :probability,
                                xlabel = "Value for $curr_param",
                                ylabel = "Probability",
                                alpha = 0.7,
                                color = :blue,
                                legend = :topright)
                if !isnothing(actual_non_time_varying_vals)
                    actual_param_value = round(getfield(actual_non_time_varying_vals, Symbol(curr_param)), digits=3)
                    scatter!(plt, [actual_param_value], [0.002], 
                            label = "Actual Value: $actual_param_value", 
                            color = :red,
                            markersize = 5,
                            marker = :circle)
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
            savefig(final_plot, plot_name_to_save)
        end
    else
        println("NO NON-TIME VARYING PARAMETER PLOTS TO DISPLAY!!!")
    end

end




function non_time_varying_param_vis(
    build_params::uciwweihr_model_params,
    data_hosp, 
    data_wastewater,
    obstimes,
    param_change_times,
    seed,
    forecast,
    forecast_weeks;
    gq_samples=nothing,
    desired_params=nothing,
    bayes_dist_type="Posterior",
    actual_non_time_varying_vals::uciwweihr_sim_params = nothing,
    save_plots::Bool=false,
    plot_name_to_save = "mcmc_nontime_varying_parameter_plots"
    )
    println("Generating non-time varying parameter plots (with priors and with wastewater)...")
    samples = uciwweihr_fit(
        data_hosp,
        data_wastewater;
        obstimes = obstimes,
        param_change_times = param_change_times,
        priors_only = true,
        seed = seed,
        params = build_params
        )
    prior_model_output = uciwweihr_gq_pp(
        samples,
        data_hosp,
        data_wastewater;
        obstimes,
        param_change_times,
        seed = seed,
        params = build_params,
        forecast = forecast, 
        forecast_weeks = forecast_weeks
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
                                title = "$bayes_dist_type $curr_param",
                                bins = 50,
                                normalize = :probability,
                                xlabel = "Value for $curr_param",
                                ylabel = "Probability",
                                alpha = 0.7,
                                color = :blue,
                                legend = :topright)
                if !isnothing(actual_non_time_varying_vals)
                    actual_param_value = round(getfield(actual_non_time_varying_vals, Symbol(curr_param)), digits=3)
                    scatter!(plt, [actual_param_value], [0.002], 
                            label = "Actual Value: $actual_param_value", 
                            color = :red,
                            markersize = 5,
                            marker = :circle)
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
            savefig(final_plot, plot_name_to_save)
        end
    else
        println("NO NON-TIME VARYING PARAMETER PLOTS TO DISPLAY!!!")
    end

end




function non_time_varying_param_vis(;
    gq_samples=nothing,
    desired_params=nothing,
    bayes_dist_type="Posterior",
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
                                title = "$bayes_dist_type $curr_param",
                                bins = 50,
                                normalize = :probability,
                                xlabel = "Probability",
                                ylabel = "Value for $curr_param",
                                xguidefont = font(8),
                                yguidefont = font(8), 
                                titlefont = font(10), 
                                legendfont = font(8),
                                legend = :topright)
                if !isnothing(actual_non_time_varying_vals)
                    actual_param_value = round(getfield(actual_non_time_varying_vals, Symbol(curr_param)), digits=3)
                    scatter!(plt, Float64[actual_param_value], [0.002], 
                            label = "Actual Value : $actual_param_value", 
                            color = :red,
                            markersize = 5,
                            marker = :circle,
                            legendfont = font(8))
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
