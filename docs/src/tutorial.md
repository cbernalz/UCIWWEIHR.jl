```@setup tutorial
using Plots; gr()
Plots.reset_defaults()
```
# [Tutorials](@id tutorial)

Welcome to the Tutorials page for the UCIWWEIHR.jl project. This section provides step-by-step guides and examples to help you get started with using the package and understanding its features.

## 1. Getting Started

### 1.1 Installation

To install the UCIWWEIHR.jl package, open the Julia REPL and run the following command:

``` julia
using Pkg
Pkg.add("UCIWWEIHR")
```

## 2. Generating simulated data

This package provides a way to also simulate data using the model specified in the future paper.  The function called `generate_simulation_data_ww_eihr` can be used to generate synthetic data for a given number of samples and features.  Here we provide a demonstration using the default settings of `generate_simulation_data_ww_eihr`:

``` @example tutorial
using UCIWWEIHR
using Plots
# Running simulation function with defaults
df = generate_simulation_data_ww_eihr()
first(df, 5)
```

Here we can make simple plots to visualize the data generated using the [Plots.jl](https://github.com/JuliaPlots/Plots.jl) package.

### 2.1. Concentration of Pathogen Genome in WW
```@example tutorial
plot(df.obstimes, df.log_ww_conc,
    label=nothing,
    xlabel="Obstimes", 
    ylabel="Conc. of Pathogen Genome in WW", 
    title="Plot of Conc. of Pathogen Genome in WW Over Time")
```

### 2.2. Hospitalizations
```@example tutorial
plot(df.obstimes, df.hosp, 
    label=nothing,
    xlabel="Obstimes", 
    ylabel="Hosp", 
    title="Plot of Hosp Over Time")
```

### 2.3. Reproductive Number
```@example tutorial
plot(df.obstimes, df.rt, 
    label=nothing,
    xlabel="Obstimes", 
    ylabel="RT", 
    title="Plot of RT Over Time")
```
