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