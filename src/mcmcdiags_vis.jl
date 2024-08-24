"""
    mcmcdiags_vis(...)
Default visualizer for results of the UCIWWEIHR model, includes posterior/priors of generated quantities and posterior predictive samples for forecasting.  Forecasting plots will have the observed data alongside.

# Arguments
- `gq_samples`: Generated quantities samples from the posterior/prior distribution, index 2 in uciwweihr_gq_pp output.
- `desired_params`: A list of lists of parameters to visualize. Each list will be visualized in a separate plot. Default is [["E_init", "I_init", "H_init"], ["gamma", "nu", "epsilon"], ["rho_gene", "tau", "df"], ["sigma_hosp"]].
- `actual_non_time_varying_vals::uciwweihr_sim_params`: A uciwweihr_sim_params object of actual non-time varying parameter values if user has access to them. Default is nothing.
- `save_plots::Bool=false`: A boolean to indicate if user wants to save the plots as pngs into a plots folder.
"""
function mcmcdiags_vis(;
    gq_samples=nothing,
    desired_params=[
        ["E_init", "I_init", "H_init"],
        ["gamma", "nu", "epsilon"],
        ["rt_init", "w_init"],
        ["rho_gene", "tau", "df"],
        ["sigma_hosp"]
    ],
    actual_non_time_varying_vals::uciwweihr_sim_params = nothing,
    save_plots::Bool=false
    )

    # Posterior/Prior Samples
    ## MCMC evaluation
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
                    plt = plot(df_param.iteration, df_param.value, 
                        label = "Chain $chain",
                        title = "Trace of $param", 
                        xlabel = "Iteration", ylabel="Value Drawn",
                        color = :black, lw = 2,
                        xguidefont = font(8),
                        yguidefont = font(8), 
                        titlefont = font(10), 
                        legendfont = font(8),
                        legend = :topright
                    )            
                    push!(cat_plots, plt) 

                    if !isnothing(actual_non_time_varying_vals)
                        actual_param_value = round(getfield(actual_non_time_varying_vals, Symbol(param)), digits=3)
                        scatter!(plt, [1], Float64[actual_param_value],
                                label = "Actual Value : $actual_param_value", 
                                color = :red,
                                markersize = 5,
                                marker = :circle,
                                legendfont = font(8))
                    end
                else
                    println("PARAMETER $param NOT IN GENERATED QUANTITIES!!!")
                end
            end
        end
    end
    if !isempty(cat_plots)
        num_plots = length(cat_plots)
        num_chains = length(unique(gq_samples.chain))
        num_params = length(desired_params)
        layout_rows = num_chains * length(desired_params[1])
        layout_cols = num_params
        
        final_plot = plot(cat_plots..., 
                            layout = (layout_rows, layout_cols), 
                            size = (1500, 1500))
        display(final_plot)
        if save_plots
            save_plots_to_docs(final_plot, "mcmc_diagnosis_plots")
        end
    else
        println("NO MCMC DIAGNOSIS PLOTS TO DISPLAY!!!")
    end

end
