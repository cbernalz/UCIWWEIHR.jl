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
        cat_plots = []
        for chain in unique(gq_samples.chain)
            for param_group in desired_params
                eff_sample_sizes = Dict{String,  Float64}()
                for param in param_group
                    if param in names(gq_samples)
                    size_temp = round(ess(gq_samples[gq_samples.chain .== chain, param]))
                    eff_sample_sizes[param] = size_temp
                    println("Effective Sample Size for $param for Chain $chain: $size_temp") 
                    else
                        println("PARAMETER $param NOT IN GENERATED QUANTITIES!!!")
                    end
                end
                long_df = stack(gq_samples[gq_samples.chain .== chain, :], Not([:iteration, :chain]), variable_name=:name, value_name=:value)
                df_filtered = filter(row -> row.name in param_group, long_df)

                for param in param_group
                    if param in names(gq_samples)
                        df_param = filter(row -> row.name == param, df_filtered)
                        title = "MCMC Diagnosis Plot for Chain $chain, $param (ESS: $(eff_sample_sizes[param]))"
                        plt = plot(df_param.iteration, df_param.value, 
                            legend=false,
                            title=title, 
                            xlabel="Iteration", ylabel="Value Drawn",
                            color = :black, lw = 2
                        )            
                        push!(cat_plots, plt) 
                    else
                        println("PARAMETER $param NOT IN GENERATED QUANTITIES!!!")
                    end
                end
            end
        end
        if !isempty(cat_plots)
            plt = plot(cat_plots..., layout=(length(unique(gq_samples.chain)) * length(desired_params[1]), length(desired_params)), size = (1000, 1000))
            display(plt)
            if save_plots
                save_plots_to_docs(plt, "mcmc_diagnosis_plots")
            end
        else
            println("NO PLOTS TO DISPLAY!!!")
        end
    else
        println("MCMC Diagnostics Plots are not requested.")
    end

    if time_varying_plots

        # Plotting time varying parameters
        var_prefixs = time_varying_params
        time_varying_plots = []
        column_names = names(gq_samples)
        for var_prefix in var_prefixs
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
            plt = plot(time_varying_plots..., layout = (length(var_prefixs), length(chains)), size = (1000, 1000))
            display(plt)
            if save_plots
                save_plots_to_docs(plt, "mcmc_time_varying_parameter_plots")
            end
        else
            println("NO TIME VARYING PARAMETER PLOTS TO DISPLAY!!!")
        end
    else
        println("MCMC time varying parameter results are not requested.")
    end



end
