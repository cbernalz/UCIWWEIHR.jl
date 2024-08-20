# [Generating simulated data with UCIWWEIHR ODE compartmental based model.](@id uciwweihr_simulation_data)

This package provides a way to also simulate data using the UCIWWEIHR ODE compartmental based model specified in the future paper.  The function called `generate_simulation_data_uciwweihr.jl` can be used to generate synthetic data for a given number of samples and features.  Here we provide a demonstration using the default settings of `generate_simulation_data_uciwweihr.jl` :

## 1. Functionality.

``` @example tutorial
using UCIWWEIHR
using Plots
# Running simulation function with defaults
params = create_uciwweihr_params()
df = generate_simulation_data_uciwweihr(params)
first(df, 5)
```

## 1.2 Visualizing UCIWWEIHR model results.

Here we can make simple plots to visualize the data generated using the [Plots](https://docs.juliaplots.org/stable/) package.

### 1.2.1. Concentration of pathogen genome in wastewater(WW).
```@example tutorial
plot(df.obstimes, df.log_ww_conc,
    label=nothing,
    xlabel="Obstimes", 
    ylabel="Conc. of Pathogen Genome in WW", 
    title="Plot of Conc. of Pathogen Genome in WW Over Time")
```

### 1.2.2. Hospitalizations.
```@example tutorial
plot(df.obstimes, df.hosp, 
    label=nothing,
    xlabel="Obstimes", 
    ylabel="Hosp", 
    title="Plot of Hosp Over Time")
```

### 1.2.3. Reproductive number.
```@example tutorial
plot(df.obstimes, df.rt, 
    label=nothing,
    xlabel="Obstimes", 
    ylabel="Rt", 
    title="Plot of Rt Over Time")
```

### 1.2.4. Hospitalization rate.
```@example tutorial
plot(df.obstimes, df.wt, 
    label=nothing,
    xlabel="Obstimes", 
    ylabel="Rt", 
    title="Plot of Hospitalization Rate Over Time")
```

## 2. Alternate Functionality.
We can also use a prespecified effective repordcution number curve or a prespecified hospitaliation probability curve.  Any combintation of presepcified or random walk curves can be used.  Here we provide an example of using both a prespecified effective reproduction number curve and a prespecified hospitalization probability curve.

``` @example tutorial
using UCIWWEIHR
using Plots
# Running simulation function with prespecified Rt and hospitalization probability
rt_custom = vcat(
    range(1, stop=1.8, length=7*4),
    fill(1.8, 7*2),
    range(1.8, stop=1, length=7*8),
    range(0.98, stop=0.8, length=7*2),
    range(0.8, stop=1.1, length=7*6),
    range(1.1, stop=0.97, length=7*3)
)
w_custom = vcat(
    range(0.3, stop=0.38, length=7*5),
    fill(0.38, 7*2),
    range(0.38, stop=0.25, length=7*8),
    range(0.25, stop=0.28, length=7*2),
    range(0.28, stop=0.34, length=7*6),
    range(0.34, stop=0.28, length=7*2)
)
params = create_uciwweihr_params(
    time_points = length(rt_custom),
    Rt = rt_custom, 
    w = w_custom
)
df = generate_simulation_data_uciwweihr(params)
first(df, 5)
```

## 2.2 Visualizing UCIWWEIHR model results.

We can visualize these results using the [Plots](https://docs.juliaplots.org/stable/) package.

### 2.2.1. Concentration of pathogen genome in wastewater(WW).
```@example tutorial
plot(df.obstimes, df.log_ww_conc,
    label=nothing,
    xlabel="Obstimes", 
    ylabel="Conc. of Pathogen Genome in WW", 
    title="Plot of Conc. of Pathogen Genome in WW Over Time")
```

### 2.2.2. Hospitalizations.
```@example tutorial
plot(df.obstimes, df.hosp, 
    label=nothing,
    xlabel="Obstimes", 
    ylabel="Hosp", 
    title="Plot of Hosp Over Time")
```

### 2.2.3. Reproductive number.
```@example tutorial
plot(df.obstimes, df.rt, 
    label=nothing,
    xlabel="Obstimes", 
    ylabel="Rt", 
    title="Plot of Rt Over Time")
```

### 2.2.4. Hospitalization rate.
```@example tutorial
plot(df.obstimes, df.wt, 
    label=nothing,
    xlabel="Obstimes", 
    ylabel="Rt", 
    title="Plot of Hospitalization Rate Over Time")
```


### [Tutorial Contents](@ref tutorial_home)