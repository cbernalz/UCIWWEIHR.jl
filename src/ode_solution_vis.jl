"""
    ode_solution_vis(...)
Visualizer for ODE solutions produced by the UCIWWEIHR model.  If actual ode solutions are provided, they will be included in the plots.

# Arguments
- `gq_samples`: Generated quantities samples from the posterior/prior distribution, index 2 in uciwweihr_gq_pp output.
- `actual_non_time_varying_vals::uciwweihr_sim_params`: A uciwweihr_sim_params object of actual non-time varying parameter values if user has access to them. Default is nothing.
- `save_plots::Bool=false`: A boolean to indicate if user wants to save the plots as pngs into a plots folder.
- `plot_name_to_save`: A string to indicate the name of the plot to save. Default is "mcmc_diagnosis_plots".
"""
function ode_solution_vis(;
    gq_samples=nothing,
    actual_E_ode_sol = nothing,
    actual_I_ode_sol = nothing,
    actual_H_ode_sol = nothing,
    quantiles = [0.5, 0.8, 0.95],
    save_plots::Bool=false,
    plot_name_to_save = "mcmc_diagnosis_plots"
    )

    ode_names = [r"E(\[|\])", r"I(\[|\])", r"H(\[|\])"]
    var_prefixes = ["E", "I", "H"]
    plots = []

    for curr_ode_name in ode_names
        var_prefix = var_prefixes[findfirst(x -> x == curr_ode_name, ode_names)]
        
        sol_df = subset_desired_ode_from_gq(gq_samples, curr_ode_name)
        medians, lower_bounds, upper_bounds = calculate_quantiles_without_chain(sol_df, var_prefix, quantiles)
        ribbon_colors = generate_colors(length(quantiles))
        time_index = 1:length(medians)
        plt = plot(title = "ODE Quantilized Solutions for $var_prefix",
                xlabel = "Time Points (daily scale)",
                ylabel = "Value for $var_prefix")
        for (i, q) in enumerate(quantiles)
            preped_upper_bounds_temp = map(x -> x[i], upper_bounds)
            preped_lower_bounds_temp = map(x -> x[i], lower_bounds)
            plot!(plt, time_index, medians, ribbon = (preped_upper_bounds_temp .- medians, medians .- preped_lower_bounds_temp),
                fillalpha = 0.2,
                label = "$(@sprintf("%.0f", q*100))% Quantile",
                color = ribbon_colors[i],
                fillcolor = ribbon_colors[i])
        end
        plot!(plt, time_index, medians, label = "Median", color = :black, lw = 2)

        if !isnothing(actual_E_ode_sol) && var_prefix == "E"
            plot!(plt, 1:length(actual_E_ode_sol), actual_E_ode_sol, label = "Actual E ODE Values", color = :red, lw = 3)
        end
        if !isnothing(actual_I_ode_sol) && var_prefix == "I"
            plot!(plt, 1:length(actual_I_ode_sol), actual_I_ode_sol, label = "Actual I ODE Values", color = :red, lw = 3)
        end
        if !isnothing(actual_H_ode_sol) && var_prefix == "H"
            plot!(plt, 1:length(actual_H_ode_sol), actual_H_ode_sol, label = "Actual H ODE Values", color = :red, lw = 3)
        end

        push!(plots, plt) 
    end

    
    if !isempty(plots)
        plt = plot(plots..., layout = (length(plots), 1), size = (1000, 1000))
        display(plt)
        if save_plots
            save_plots_to_docs(plt, plot_name_to_save)
        end
    else
        println("NO ODE SOLUTION PLOTS TO DISPLAY!!!")
    end


end
