"""
    optimize_many_MAP(model, n_reps = 100, top_n = 1, verbose = true)

Try n_reps different initializations to get MAP estimate.

Function by Damon Bayer

"""
function optimize_many_MAP(model, n_reps = 100, top_n = 1, verbose = true)
  lp_res = repeat([-Inf], n_reps)
  for i in eachindex(lp_res)
      if verbose
          println(i)
      end
      Random.seed!(i)
      try
          lp_res[i] = optimize(model, MAP(), LBFGS(linesearch = LineSearches.BackTracking())).lp
      catch
      end
  end
  eligible_indices = findall(.!isnan.(lp_res) .& isfinite.(lp_res))
  best_n_seeds =  eligible_indices[sortperm(lp_res[eligible_indices], rev = true)][1:top_n]

  map(best_n_seeds) do seed
    Random.seed!(seed)
    optimize(model, MAP(), LBFGS(linesearch = LineSearches.BackTracking())).values.array
  end
end


"""
    optimize_many_MAP2(model, n_reps = 100, top_n = 1, verbose = true)

Try n_reps different initializations to get MAP estimate.

Modified by Christian Bernal Zelaya
"""
function optimize_many_MAP2(model, n_reps=100, top_n=1, verbose=true)
    println("Optimizing initializations....")
    lp_res = repeat([-Inf], n_reps)
    for i in eachindex(lp_res)
        if verbose
            println("Trial: $i")
        end
        Random.seed!(i)
        try
            result = optimize(model, MAP(), LBFGS(linesearch=LineSearches.BackTracking()))
            lp_res[i] = result.lp
            if verbose
                println("Optimization successful for trial $i with log-probability: $(result.lp)")
            end
        catch e
            if verbose
                println("Optimization failed for trial $i")
            end
        end
    end
    eligible_indices = findall(.!isnan.(lp_res) .& isfinite.(lp_res))
    best_n_seeds = eligible_indices[sortperm(lp_res[eligible_indices], rev=true)][1:top_n]
    map(best_n_seeds) do seed
        Random.seed!(seed)
        result = optimize(model, MAP(), LBFGS(linesearch=LineSearches.BackTracking()))
        println("Optimization found for seed $seed with log-probability: $(result.lp)")
        return result.values.array
    end
end

"""
    optimize_many_MAP2_wrapper(...)

Wrapper function for optimize_many_MAP2 that uses model based on inputs to function.

Created by Christian Bernal Zelaya
"""
function optimize_many_MAP2_wrapper(
    data_hosp,
    data_wastewater,
    obstimes_hosp,
    obstimes_wastewater,
    param_change_times,
    params::model_params_time_var_hosp_prev;
    n_reps=100,
    top_n=1,
    verbose=true,
    warning_bool=true
)
# prevalence model
    println("Optimizing initializations for uciwweihr_model with wastewater and time-varying hospitalization probability - Prevalence Model!!!")
    ## model with wastewater and time-varying hospitalization probability
    obstimes_hosp = convert(Vector{Int64}, obstimes_hosp)
    obstimes_wastewater = convert(Vector{Int64}, obstimes_wastewater)
    param_change_times = convert(Vector{Int64}, param_change_times)
    my_model = uciwweihr_model(
        data_hosp,
        data_wastewater,
        obstimes_hosp,
        obstimes_wastewater,
        param_change_times,
        params;
        warning_bool=warning_bool
    )
    return optimize_many_MAP2(my_model, n_reps, top_n, verbose)
end


function optimize_many_MAP2_wrapper(
    data_hosp,
    data_wastewater,
    obstimes_hosp,
    obstimes_wastewater,
    param_change_times,
    params::model_params_time_var_hosp_inc;
    n_reps=100,
    top_n=1,
    verbose=true,
    warning_bool=true
)
# incidence model
    println("Optimizing initializations for uciwweihr_model with wastewater and time-varying hospitalization probability - Incidence Model!!!")
    ## model with wastewater and time-varying hospitalization probability
    obstimes_hosp = convert(Vector{Int64}, obstimes_hosp)
    obstimes_wastewater = convert(Vector{Int64}, obstimes_wastewater)
    param_change_times = convert(Vector{Int64}, param_change_times)
    my_model = uciwweihr_model(
        data_hosp,
        data_wastewater,
        obstimes_hosp,
        obstimes_wastewater,
        param_change_times,
        params;
        warning_bool=warning_bool
    )
    return optimize_many_MAP2(my_model, n_reps, top_n, verbose)
end

function optimize_many_MAP2_wrapper(
    data_hosp,
    obstimes_hosp,
    param_change_times,
    params::model_params_non_time_var_hosp_no_ww;
    n_reps=100,
    top_n=1,
    verbose=true,
    warning_bool=true
)

    println("Optimizing initializations for uciwweihr_model without wastewater and nontime-varying hospitalization probability!!!")
    ## model without wastewater and without time-varying hospitalization probability
    obstimes_hosp = convert(Vector{Int64}, obstimes_hosp)
    param_change_times = convert(Vector{Int64}, param_change_times)
    my_model = uciwweihr_model(
        data_hosp,
        obstimes_hosp,
        param_change_times,
        params;
        warning_bool=warning_bool
    )
    return optimize_many_MAP2(my_model, n_reps, top_n, verbose)
end


## model with wastewater and non-time-varying hospitalization probability - not implemented
## model without wastewater and with time-varying hospitalization probability - not implemented