# [Generating simulated data with UCIWWEIHR ODE compartmental based model.](@id uciwweihr_simulation_data)

This package provides a way to also simulate data using the UCIWWEIHR ODE compartmental based model specified in the future paper.  The function called `generate_simulation_data_uciwweihr.jl` can be used to generate synthetic data for a given number of samples and features.  Here we provide a demonstration using the default settings of `generate_simulation_data_uciwweihr.jl` :

## 1. Functionality.

``` @example tutorial
using UCIWWEIHR
using Plots
# Running simulation function with defaults
df = generate_simulation_data_uciwweihr()
first(df, 5)
```

## 2. Visualizing UCIWWEIHR model results.

Here we can make simple plots to visualize the data generated using the [TidierPlots](https://tidierorg.github.io/TidierPlots.jl/stable/) package.

### 2.1. Concentration of pathogen genome in wastewater(WW).
```@example tutorial
plot(df.obstimes, df.log_ww_conc,
    label=nothing,
    xlabel="Obstimes", 
    ylabel="Conc. of Pathogen Genome in WW", 
    title="Plot of Conc. of Pathogen Genome in WW Over Time")
```

### 2.2. Hospitalizations.
```@example tutorial
plot(df.obstimes, df.hosp, 
    label=nothing,
    xlabel="Obstimes", 
    ylabel="Hosp", 
    title="Plot of Hosp Over Time")
```

### 2.3. Reproductive number.
```@example tutorial
plot(df.obstimes, df.rt, 
    label=nothing,
    xlabel="Obstimes", 
    ylabel="Rt", 
    title="Plot of Rt Over Time")
```


### [Tutorial Contents](@ref tutorial_home)