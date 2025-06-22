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
    time_points = 70
)
df = generate_simulation_data_uciwweihr(params)

params_ext = create_uciwweihr_sim_params(
    time_points = 84
)
df_ext = generate_simulation_data_uciwweihr(params_ext)
first(df, 5)
```

``` @example tutorial_forecast
first(df_ext, 5)
```

## 2. Sampling from the Posterior Distribution and Posterior Predictive Distribution.

Here we sample from the posterior distribution using the `uciwweihr_fit.jl` function.  First, we setup some presets, where we need to use `create_uciwweihr_model_params()` to get default parameters for the model.  Then we have an array where index 1 contains the posterior/prior predictive samples, index 2 contains the posterior/prior generated quantities samples, and index 3 contains the original sampled parameters for the model.  The difference here is that we set `forecast = true` and `forecast_weeks = 4` to forecast 4 weeks into the future.  One other thing to note, is that we allow misalignment of hospital and wastewater data's observed times.  For this tutorial, we use the same observed points.

``` @example tutorial_forecast
data_hosp = df.hosp
data_wastewater = df.log_ww_conc
obstimes_hosp = df.obstimes
obstimes_wastewater = df.obstimes
max_obstime = max(length(obstimes_hosp), length(obstimes_wastewater))
param_change_times = 1:7:max_obstime # Change every week
priors_only = false
n_samples = 500
forecast = true
forecast_days = 14

E_init_sd=0.2; log_E_init_mean=log(200)
I_init_sd=0.2; log_I_init_mean=log(100)
H_init_sd=0.2; log_H_init_mean=log(20)
gamma_sd=0.02; log_gamma_mean=log(1/4)
nu_sd=0.02; log_nu_mean=log(1/7)
epsilon_sd=0.02; log_epsilon_mean=log(1/5)
rho_gene_sd=0.02; log_rho_gene_mean=log(0.011)
    
sigma_ww_sd=0.02; log_sigma_ww_mean=log(0.1)
sigma_hosp_sd=0.01; log_sigma_hosp_mean=log(500.0)

Rt_init_sd=0.3; Rt_init_mean=0.2
sigma_Rt_sd=0.2; sigma_Rt_mean=-3.0
w_init_sd=0.04; w_init_mean=logit(0.35)
sigma_w_sd=0.2; sigma_w_mean=-3.5
message = true
model_params = create_model_params_time_var_hosp(
    E_init_sd, log_E_init_mean,
    I_init_sd, log_I_init_mean,
    H_init_sd, log_H_init_mean,
    gamma_sd, log_gamma_mean,
    nu_sd, log_nu_mean,
    epsilon_sd, log_epsilon_mean,
    rho_gene_sd, log_rho_gene_mean,
    sigma_ww_sd, log_sigma_ww_mean,
    sigma_hosp_sd, log_sigma_hosp_mean,
    Rt_init_sd, Rt_init_mean,
    sigma_Rt_sd, sigma_Rt_mean,
    w_init_sd, w_init_mean,
    sigma_w_sd, sigma_w_mean,
    message;
)
init_params = optimize_many_MAP2_wrapper(
    data_hosp,
    data_wastewater,
    obstimes_hosp,
    obstimes_wastewater,
    param_change_times,
    model_params;
    verbose=false,
    warning_bool=false,
)
samples = fit(
    data_hosp,
    data_wastewater,
    obstimes_hosp,
    obstimes_wastewater,
    param_change_times,
    model_params;
    priors_only,
    n_samples,
    n_discard_initial = 200,
    n_chains = 1,
    init_params = init_params
    )
model_output = generate_pq_pp(
    samples,
    data_hosp,
    data_wastewater,
    obstimes_hosp,
    obstimes_wastewater,
    param_change_times,
    model_params;
    forecast=true, forecast_days=forecast_days
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

We can again look at model diagnostics, posterior distribution of time or non-time varying parameters, and the posterior predictive distribution extended for forecasting.  We can also add certain parameters to ensure priors will be plotted alongside their corresponding posteriors.

```@example tutorial_forecast
uciwweihr_visualizer(
    data_hosp, 
    data_wastewater,
    forecast_days,
    obstimes_hosp,
    obstimes_wastewater,
    param_change_times,
    2024,
    forecast,
    model_params;
    pp_samples = model_output[1],
    gq_samples = model_output[2],
    samples = model_output[3],
    obs_data_hosp = df_ext.hosp,
    obs_data_wastewater = df_ext.log_ww_conc, 
    actual_rt_vals = df_ext.rt, 
    actual_w_t = df_ext.wt, 
    actual_E_ode_sol = df.E_ode_comp_sol,
    actual_I_ode_sol = df.I_ode_comp_sol,
    actual_H_ode_sol = df.H_ode_comp_sol,
    actual_non_time_varying_vals = params,
    bayes_dist_type = "Posterior",
    save_plots = true,
    plot_name_to_save_mcmcdiag = "plots/mcmc_diagnosis_plots1",
    plot_name_to_save_time_varying = "plots/mcmc_time_varying_parameter_plots1",
    plot_name_to_save_non_time_varying = "plots/mcmc_nontime_varying_parameter_plots1",
    plot_name_to_save_ode_sol = "plots/mcmc_ode_solution_plots1",
    plot_name_to_save_pred_param = "plots/mcmc_pred_parameter_plots1",
    plot_name_to_save_log_like = "plots/mcmc_log_prob_trace_plot1"
)
```

### 3.1. MCMC Diagnostic Plots.

![Plot 1](plots/mcmc_diagnosis_plots1.png)

### 3.2. Time Varying Parameter Results Plot.

![Plot 2](plots/mcmc_time_varying_parameter_plots1.png)

### 3.3. Non-Time Varying Parameter Results Plot.
![Plot 3](plots/mcmc_nontime_varying_parameter_plots1.png)

### 3.4. ODE Solution Plot.
![Plot 4](plots/mcmc_ode_solution_plots1.png)

### 3.5. Posterior Predictive Distribution Plot.

![Plot 4](plots/mcmc_pred_parameter_plots1.png)

### 3.6. Log Prob Trace Plot.
![Plot 5](plots/mcmc_log_prob_trace_plot1.png)


### [Tutorial Contents](@ref tutorial_home)