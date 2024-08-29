```@setup tutorial_forecast
using Plots, StatsPlots; gr()
Plots.reset_defaults()

```

# [Generating Posterior Distribution Samples with UCIWWEIHR ODE Compartmental Based Model with Forecasting.](@id uciwwiehr_model_fitting_with_forecast)

Here we extend the [previous tutorial](@ref uciwwiehr_model_fitting_no_forecast)  to include forecasting capabilities.  We start with generating out data using `generate_simulation_data_uciwweihr`'s alternate parameterization where we do not prespecify the effective reproduction number and hospitalization probability but instead preform a log-normal random walk and a logit-normal random walk respectively.  We then sample from the posterior distribution using the `uciwweihr_fit.jl` function.  We then generate desired quantities and forecast for a given time period with the posterior predictive distribution, using `uciwweihr_gq_pp.jl`.


## 1. Data Generation.

Here we generate two datasets, one with 150 time points and one with 178 time points.  We will use the 150 time point dataset for fitting and the 178 time point dataset for forecast evaluation.

``` @example tutorial_forecast
using UCIWWEIHR
# Running simulation function with presets
params = create_uciwweihr_sim_params(
    time_points = 150
)
df = generate_simulation_data_uciwweihr(params)

params_ext = create_uciwweihr_sim_params(
    time_points = 178
)
df_ext = generate_simulation_data_uciwweihr(params_ext)
first(df, 5)
```

``` @example tutorial_forecast
first(df_ext, 5)
```

## 2. Sampling from the Posterior Distribution and Posterior Predictive Distribution.

Here we sample from the posterior distribution using the `uciwweihr_fit.jl` function.  First, we setup some presets, where we need to use `create_uciwweihr_model_params()` to get default parameters for the model.  Then we have an array where index 1 contains the posterior/prior predictive samples, index 2 contains the posterior/prior generated quantities samples, and index 3 contains the original sampled parameters for the model.  The difference here is that we set `forecast = true` and `forecast_weeks = 4` to forecast 4 weeks into the future.

``` @example tutorial_forecast
data_hosp = df.hosp
data_wastewater = df.log_ww_conc
obstimes = df.obstimes
param_change_times = 1:7:length(obstimes) # Change every week
priors_only = false
n_samples = 200
forecast = true
forecast_weeks = 4

model_params = create_uciwweihr_model_params()
samples = uciwweihr_fit(
    data_hosp,
    data_wastewater;
    obstimes,
    param_change_times,
    priors_only,
    n_samples,
    params = model_params
)
model_output = uciwweihr_gq_pp(
    samples,
    data_hosp,
    data_wastewater;
    obstimes = obstimes,
    param_change_times = param_change_times,
    params = model_params,
    forecast = forecast,
    forecast_weeks = forecast_weeks
)

first(model_output[1][:,1:5], 5)
```

``` @example tutorial_forecast
first(model_output[2][:,1:5], 5)
```

``` @example tutorial_forecast
first(model_output[3][:,1:5], 5)
```

## 3. MCMC Diagnostic Plots/Results Along with Posterior Predictive Distribution.

We can again look at model diagnostics, posterior distribution of time or non-time varying parameters, and the posterior predictive distribution extended for forecasting.

```@example tutorial_forecast
uciwweihr_visualizer(
    pp_samples = model_output[1],
    gq_samples = model_output[2],
    data_hosp = df_ext.hosp,
    data_wastewater = df_ext.log_ww_conc, 
    actual_rt_vals = df_ext.rt, 
    actual_w_t = df_ext.wt, 
    actual_non_time_varying_vals = params,
    forecast_weeks = forecast_weeks,
    bayes_dist_type = "Posterior",
    save_plots = true,
    plot_name_to_save_mcmcdiag = "mcmc_diagnosis_plots1",
    plot_name_to_save_time_varying = "mcmc_time_varying_parameter_plots1",
    plot_name_to_save_non_time_varying = "mcmc_nontime_varying_parameter_plots1",
    plot_name_to_save_pred_param = "mcmc_pred_parameter_plots1"
)
```

### 3.1. MCMC Diagnostic Plots.

![Plot 1](plots/mcmc_diagnosis_plots1.png)

### 3.2. Time Varying Parameter Results Plot.

![Plot 2](plots/mcmc_time_varying_parameter_plots1.png)

### 3.3. Non-Time Varying Parameter Results Plot.
![Plot 3](plots/mcmc_nontime_varying_parameter_plots1.png)

### 3.4. Posterior Predictive Distribution Plot.

![Plot 4](plots/mcmc_pred_parameter_plots1.png)


### [Tutorial Contents](@ref tutorial_home)