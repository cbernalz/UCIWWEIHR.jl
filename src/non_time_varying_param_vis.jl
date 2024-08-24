"""
    non_time_varying_param_vis(...)

Used in the `uciwweihr_visualizer` to create visuals for non-time varying parameters.

# Arguments
- `gq_samples`: Generated quantities samples from the posterior/prior distribution, index 2 in uciwweihr_gq_pp output.
- `desired_params`: A list of lists of parameters to visualize. Each list will be visualized in a separate plot. Default is any parameter not in this list : ["alpha_t", "w_t", "rt_vals", "log_genes_mean", "H"]
- `bayes_dist_type`: A string to indicate if user is using Posterior or Prior distribution. Default is "Posterior".
- `actual_non_time_varying_vals::uciwweihr_sim_params`: A uciwweihr_sim_params object of actual non-time varying parameter values if user has access to them. Default is nothing.
- `save_plots::Bool=false`: A boolean to indicate if user wants to save the plots as pngs into a plots folder.
"""
function non_time_varying_param_vis(;
    gq_samples=nothing,
    desired_params=nothing,
    bayes_dist_type="Posterior",
    actual_non_time_varying_vals::uciwweihr_sim_params = nothing,
    save_plots::Bool=false
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
            save_plots_to_docs(final_plot, "mcmc_nontime_varying_parameter_plots")
        end
    else
        println("NO NON-TIME VARYING PARAMETER PLOTS TO DISPLAY!!!")
    end
    
    
end
