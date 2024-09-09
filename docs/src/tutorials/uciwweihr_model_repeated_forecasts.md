```@setup tutorial_forecast
using Plots, StatsPlots; gr()
Plots.reset_defaults()

```

# [Generating Repeated Forecasts Using the UCIWWEIHR model.](@id uciwwiehr_model_repeated_forecasts)

Here we show how we can construct repeated forecasts using the UCIWWEIHR model.  We start with generating out data using `generate_simulation_data_uciwweihr`'s alternate parameterization where we do prespecify the effective reproduction number and hospitalization probability.  



## 1. Data Generation.

Here we simulate a dataset, one with 175 time points.  

``` @example tutorial_forecast
using UCIWWEIHR
# Running simulation function with presets
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
params = create_uciwweihr_sim_params(
    time_points = length(rt_custom),
    Rt = rt_custom, 
    w = w_custom
)
df = generate_simulation_data_uciwweihr(params)
```

## 2. Constructing Repeat Forecasts.

We use the `repeated_forecast` function to generate forecasts for a given number of weeks, for a given number of time points.  Along with this we need to specify presets.  Output of this function is an array with the first index controlling which result we are looking at.  The next contains a `uciwweihr_gq_pp` output.

``` @example tutorial_forecast
data_hosp = df.hosp
data_wastewater = df.log_ww_conc
obstimes_hosp = df.obstimes
obstimes_wastewater = df.obstimes
max_obstime = max(length(obstimes_hosp), length(obstimes_wastewater))
param_change_times = 1:7:max_obstime # Change every week
priors_only = false
n_samples = 200
n_forecast_weeks = 2
forecast_points = [
    param_change_times[end-5],
    param_change_times[end-4],
    param_change_times[end-3],
    param_change_times[end-2]
]

model_params = create_uciwweihr_model_params()

rep_results = repeated_forecast(
    data_hosp,
    data_wastewater,
    obstimes_hosp,
    obstimes_wastewater;
    n_samples = n_samples,
    params = model_params,
    n_forecast_weeks = 2,
    forecast_points = forecast_points
)
```

## 3. Visualizing Results Of Repeated Forecasts.

We can take a look at these forecasts using the `uciwweihr_visualizer` function.  We can also add certain parameters to ensure we only see the plots we want.

```@example tutorial_forecast
for res_index in 1:length(forecast_points)
    uciwweihr_visualizer(
        data_hosp, 
        data_wastewater,
        n_forecast_weeks,
        obstimes_hosp,
        obstimes_wastewater,
        param_change_times,
        2024,
        true,
        model_params;
        pp_samples = rep_results[res_index][2][1],
        gq_samples = rep_results[res_index][2][2],
        obs_data_hosp = data_hosp,
        obs_data_wastewater = data_wastewater, 
        actual_rt_vals = df.rt, 
        actual_w_t = df.wt, 
        actual_non_time_varying_vals = params,
        bayes_dist_type = "Posterior",
        mcmcdaigs = false,
        time_varying_plots = false,
        non_time_varying_plots = false,
        pred_param_plots = true,
        save_plots = true,
        plot_name_to_save_pred_param = "mcmc_pred_parameter_plots_rep_res"*string(res_index)*".png"
    )
end
```

### 3.1. Forecast Point 1.

![Plot 1](plots/mcmc_pred_parameter_plots_rep_res1.png)

### 3.2. Forecast Point 2.

![Plot 2](plots/mcmc_pred_parameter_plots_rep_res2.png)

### 3.3. Forecast Point 3.

![Plot 3](plots/mcmc_pred_parameter_plots_rep_res3.png)

### 3.4. Forecast Point 4.

![Plot 4](plots/mcmc_pred_parameter_plots_rep_res4.png)


### [Tutorial Contents](@ref tutorial_home)