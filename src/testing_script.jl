using AxisArrays
using MCMCChains
using Optim
using LineSearches
using Random
using LinearAlgebra
using NaNMath
using LogExpFunctions
using PreallocationTools
using Distributions
using Turing
using Random
using OrdinaryDiffEq
using Logging
using CSV
using DataFrames
using StatsBase
using Plots
using Printf
using Colors

include("model_params.jl")
include("likelihood_helpers.jl")
include("generate_simulation_data_uciwweihr.jl")
include("fit.jl")
include("generate_pq_pp.jl")
include("eihr_ode.jl")
include("negativebinomial2.jl")
include("uciwweihr_model.jl")
include("helper_functions.jl")

# Running simulation function with presets
params = create_uciwweihr_sim_params(
    time_points = 70
)
df = generate_simulation_data_uciwweihr(params)

# params_ext = create_uciwweihr_sim_params(
#     time_points = 84
# )
# df_ext = generate_simulation_data_uciwweihr(params_ext)
first(df, 5)
# first(df_ext, 5)
data_hosp_prev = df.hosp
data_hosp_inc = df.hosp_inc
data_wastewater = df.log_ww_conc
obstimes_hosp_prev = df.obstimes
obstimes_hosp_inc = df.obstimes
obstimes_wastewater = df.obstimes
max_obstime = max(length(obstimes_hosp_prev), length(obstimes_hosp_inc), length(obstimes_wastewater))
param_change_times = 1:7:max_obstime # Change every week
priors_only = false
n_samples = 20
forecast = true
forecast_days = 14

E_init_sd=0.2; log_E_init_mean=log(200)
I_init_sd=0.2; log_I_init_mean=log(100)
H_init_sd=0.2; log_H_init_mean=log(20)
CH_init_sd=0.2; log_CH_init_mean=log(20)
gamma_sd=0.02; log_gamma_mean=log(1/4)
nu_sd=0.02; log_nu_mean=log(1/7)
epsilon_sd=0.02; log_epsilon_mean=log(1/5)
rho_gene_sd=0.02; log_rho_gene_mean=log(0.011)
    
sigma_ww_sd=0.02; log_sigma_ww_mean=log(0.1)
sigma_hosp_prev_sd=0.01; log_sigma_hosp_prev_mean=log(500.0)
sigma_hosp_inc_sd=0.01; log_sigma_hosp_inc_mean=log(500.0)

Rt_init_sd=0.3; Rt_init_mean=0.2
sigma_Rt_sd=0.2; sigma_Rt_mean=-3.0
w_init_sd=0.04; w_init_mean=logit(0.35)
sigma_w_sd=0.2; sigma_w_mean=-3.5
message = true
model_params = create_model_params_hosp_inc_prev(
    E_init_sd, log_E_init_mean,
    I_init_sd, log_I_init_mean,
    H_init_sd, log_H_init_mean,
    CH_init_sd, log_CH_init_mean,
    gamma_sd, log_gamma_mean,
    nu_sd, log_nu_mean,
    rho_gene_sd, log_rho_gene_mean,
    epsilon_sd, log_epsilon_mean,
    sigma_ww_sd, log_sigma_ww_mean,
    sigma_hosp_prev_sd, log_sigma_hosp_prev_mean,
    sigma_hosp_inc_sd, log_sigma_hosp_inc_mean,
    Rt_init_sd, Rt_init_mean,
    sigma_Rt_sd, sigma_Rt_mean,
    w_init_sd, w_init_mean,
    sigma_w_sd, sigma_w_mean,
    message;
)
init_params = [
    0.0, # E_init_non_centered
    0.0, # I_init_non_centered
    0.0, # H_init_non_centered
    0.0, # CH_init_non_centered
    0.0, # gamma_non_centered
    0.0, # nu_non_centered
    0.0, # epsilon_non_centered
    0.0, # rho_gene_non_centered
    0.0, # sigma_ww_non_centered
    0.0, # sigma_hosp_inc_non_centered
    0.0, # sigma_hosp_prev_non_centered


    # Rt_params_non_centered (11)
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0,

    # w_params_non_centered (11)
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0,
]
#TODO: Ask about likelihood_helpers needing CH_init_non_centerd (CH_init = H_init line 239)
samples = fit(
    data_hosp_prev,
    data_hosp_inc,
    data_wastewater,
    obstimes_hosp_prev,
    obstimes_hosp_inc,
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
    data_hosp_prev,
    data_hosp_inc,
    data_wastewater,
    obstimes_hosp_prev,
    obstimes_hosp_inc,
    obstimes_wastewater,
    param_change_times,
    model_params;
    forecast=true, forecast_days=forecast_days
)

first(model_output[1][:,1:5], 5)
first(model_output[2][:,1:5], 5)
first(model_output[3][:,1:5], 5)