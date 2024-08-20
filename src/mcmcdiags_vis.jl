"""
    mcmcdiags_vis(...)
Default visualizer for results of the UCIWWEIHR model, includes posterior/priors of generated quantities and posterior predictive samples for forecasting.  Forecasting plots will have the observed data alongside.

# Arguments
- `gq_samples`: Generated quantities samples from the posterior/prior distribution, index 2 in uciwweihr_gq_pp output.
- `desired_params`: A list of lists of parameters to visualize. Each list will be visualized in a separate plot. Default is [["E_init", "I_init", "H_init"], ["gamma", "nu", "epsilon"], ["rho_gene", "tau", "df"], ["sigma_hosp"]].
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

end
