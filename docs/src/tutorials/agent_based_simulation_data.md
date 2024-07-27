# [Generating simulated data with an agent based model.](@id agent_based_simulation_data)

This package provides a way to also simulate data using the agent based model in the future paper.  The function called `generate_simulation_data_agent.jl` can be used to generate synthetic data for a given population size and features.  Here we provide a demonstration using the default settings of `generate_simulation_data_agent.jl` :


## 1. Functionality.

``` @example tutorial
using UCIWWEIHR
using Plots
# Running simulation function with defaults
df = generate_simulation_data_agent()
first(df, 5)
```

## 2. Visualizing SEIHR compartments.

We can also use the [TidierPlots](https://tidierorg.github.io/TidierPlots.jl/stable/) package to visualize the data generated.

```@example tutorial
plot(df.Time, df.S, label = "Suseptible", 
    xlabel = "Time", 
    ylabel = "Number of Individuals", 
    title = "Agent Based Model Simulation Results")
plot!(df.Time, df.E, label = "Exposed")
plot!(df.Time, df.I, label = "Infected")
plot!(df.Time, df.H, label = "Hospitalized")
plot!(df.Time, df.R, label = "Recovered")
```


### [Tutorial Contents](@ref tutorial_home)