"""
    predictive_param_vis(...)
Used in the `uciwweihr_visualizer` to create visuals for wastewater data and hospitalization data.

# Arguments
- `pp_samples`: A DataFrame of posterior or prior predictive samples.
- `data_wastewater`: An array of actual wastewater values if user has access to them assumed, using time scale of observed time points.  Default is nothing.
- `data_hosp`: An array of actual hospitalization values if user has access to them assumed, , using time scale of observed time points.  Default is nothing.
- `forecast_weeks`: An integer of the number of weeks forecasted. Default is 0.
- `vars_to_pred`: A list of variables to predict. Default is ["data_wastewater", "data_hosp"].
- `quantiles`: A list of quantiles to calculate for ploting uncertainty. Default is [0.5, 0.8, 0.95].
- `bayes_dist_type`: A string to indicate if user is using Posterior or Prior distribution. Default is "Posterior".
- `save_plots::Bool=false`: A boolean to indicate if user wants to save the plots as pngs into a plots folder.
"""
function predictive_param_vis(;
    pp_samples = nothing,
    data_wastewater = nothing,
    data_hosp = nothing,
    forecast_weeks = 0,
    vars_to_pred = ["data_wastewater", "data_hosp"],
    quantiles = [0.5, 0.8, 0.95],
    bayes_dist_type = "Posterior",
    save_plots::Bool = false
    )
    # Plotting wastewater and hosp predictive posterior / prior
    pred_plots = []
    column_names = names(pp_samples)
    for var_prefix in vars_to_pred
        pred_var_names = filter(name -> startswith_any(name, [var_prefix]), column_names)
        pred_var_df = pp_samples[:, [pred_var_names..., "iteration", "chain"]]
        chains = unique(pred_var_df.chain)
        for chain in chains
            pred_var_elem = pred_var_df[:, [pred_var_names..., "chain"]]
            medians, lower_bounds, upper_bounds = calculate_quantiles(pred_var_elem, chain, var_prefix, quantiles)
            ribbon_colors = generate_colors(length(quantiles))
            preped_medians = repeat_last_n_elements(medians, forecast_weeks, 7)
            preped_lower_bounds = repeat_last_n_elements(lower_bounds, forecast_weeks, 7)
            preped_upper_bounds = repeat_last_n_elements(upper_bounds, forecast_weeks, 7)
            time_index = 1:length(preped_medians)
            plt = plot(title = "$bayes_dist_type Quantiles for Chain $chain for $var_prefix",
                        xlabel = "Time Points (daily scale)",
                        ylabel = "Value for $var_prefix")
            for (i, q) in enumerate(quantiles)
                preped_upper_bounds_temp = map(x -> x[i], preped_upper_bounds)
                preped_lower_bounds_temp = map(x -> x[i], preped_lower_bounds)
                plot!(plt, time_index, preped_medians, ribbon = (preped_upper_bounds_temp .- preped_medians, preped_medians .- preped_lower_bounds_temp),
                        fillalpha = 0.2,
                        label = "$(@sprintf("%.0f", q*100))% Quantile",
                        color = ribbon_colors[i],
                        fillcolor = ribbon_colors[i])
            end
            plot!(plt, time_index, preped_medians, label = "Median", color = :black, lw = 2)
            if !isnothing(data_wastewater) && var_prefix == "data_wastewater"
                scatter!(plt, 1:length(data_wastewater), data_wastewater, label = "Actual Wastewater Values", color = :red, lw = 2, marker = :circle)
            end
            if !isnothing(data_hosp) && var_prefix == "data_hosp"
                scatter!(plt, 1:length(data_hosp), data_hosp, label = "Actual Hosp Values", color = :red, lw = 2, marker = :circle)
            end
            push!(pred_plots, plt) 
        end
    end

    if !isempty(pred_plots)
        chains = unique(pp_samples.chain)
        plt = plot(pred_plots..., layout = (length(pred_plots), length(chains)), size = (1000, 1000))
        display(plt)
        if save_plots
            save_plots_to_docs(plt, "mcmc_pred_parameter_plots")
        end
    else
        println("NO TIME VARYING PARAMETER PLOTS TO DISPLAY!!!")
    end
    
end
