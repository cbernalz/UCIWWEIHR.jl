"""
    uciwweihr_model_params

Struct for holding parameters used in the UCIWWEIHR ODE compartmental model.  Use `create_uciwweihr_model_params` to create an instance of this struct.

# Fields
- `E_init_sd::Float64=50.0`: Standard deviation for the initial number of exposed individuals.
- `E_init_mean::Int64=200`: Mean for the initial number of exposed individuals.
- `I_init_sd::Float64=20.0`: Standard deviation for the initial number of infected individuals.
- `I_init_mean::Int64=100`: Mean for the initial number of infected individuals.
- `H_init_sd::Float64=5.0`: Standard deviation for the initial number of hospitalized individuals.
- `H_init_mean::Int64=20`: Mean for the initial number of hospitalized individuals.
- `gamma_sd::Float64=0.02`: Standard deviation for the rate of incubation.
- `log_gamma_mean::Float64=log(1/4)`: Mean for the rate of incubation on log scale.
- `nu_sd::Float64=0.02`: Standard deviation for the rate of leaving the infected compartment.
- `log_nu_mean::Float64=log(1/7)`: Mean for the rate of leaving the infected compartment on the log scale.
- `epsilon_sd::Float64=0.02`: Standard deviation for the rate of hospitalization recovery.
- `log_epsilon_mean::Float64=log(1/5)`: Mean for the rate of hospitalization recovery on the log scale.
- `rho_gene_sd::Float64=0.02`: Standard deviation for the rho prior.
- `log_rho_gene_mean::Float64=log(0.011)`: Mean for the row prior on log scale.
- `tau_sd::Float64=0.02`: Standard deviation for the scale/variation of the log scale data.
- `log_tau_mean::Float64=log(0.1)`: Mean for the scale/variation of the log scale data on log scale itself.
- `df_shape::Float64=2.0`: Shape parameter for the gamma distribution.
- `df_scale::Float64=10.0`: Scale parameter for the gamma distribution.
- `sigma_hosp_sd::Float64=50.0`: Standard deviation for the negative binomial distribution for hospital data.
- `sigma_hosp_mean::Float64=500.0`: Mean for the negative binomial distribution for hospital data.
- `Rt_init_sd::Float64=0.3`: Standard deviation for the initial value of the time-varying reproduction number.
- `Rt_init_mean::Float64=0.2`: Mean for the initial value of the time-varying reproduction number.
- `sigma_Rt_sd::Float64=0.2`: Standard deviation for normal prior of log time-varying reproduction number standard deviation.
- `sigma_Rt_mean::Float64=-3.0`: Mean for normal prior of log time-varying reproduction number standard deviation.
- `w_init_sd::Float64=0.1`: Standard deviation for the initial value of the time-varying hospitalization rate.
- `w_init_mean::Float64=log(0.35)`: Mean for the initial value of the time-varying hospitalization rate.
- `sigma_w_sd::Float64=0.2`: Standard deviation for normal prior of log time-varying hospitalization rate standard deviation.
- `sigma_w_mean::Float64=-3.5`: Mean for normal prior of time-varying hospitalization rate standard deviation.
"""
struct uciwweihr_model_params
    E_init_sd::Float64
    E_init_mean::Int64
    I_init_sd::Float64
    I_init_mean::Int64
    H_init_sd::Float64
    H_init_mean::Int64
    gamma_sd::Float64
    log_gamma_mean::Float64
    nu_sd::Float64
    log_nu_mean::Float64
    epsilon_sd::Float64
    log_epsilon_mean::Float64
    rho_gene_sd::Float64
    log_rho_gene_mean::Float64
    tau_sd::Float64
    log_tau_mean::Float64
    df_shape::Float64
    df_scale::Float64
    sigma_hosp_sd::Float64
    sigma_hosp_mean::Float64
    Rt_init_sd::Float64
    Rt_init_mean::Float64
    sigma_Rt_sd::Float64
    sigma_Rt_mean::Float64
    w_init_sd::Float64
    w_init_mean::Float64
    sigma_w_sd::Float64
    sigma_w_mean::Float64
end

"""
    create_uciwweihr_model_params(; kwargs...)

Creates a `uciwweihr_sim_params` struct with the option to either use a predetermined `Rt` and `w` or generate them as random walks.

# Arguments
- `kwargs...`: Named arguments corresponding to the fields in `uciwweihr_sim_params`.

# Returns
- `params::uciwweihr_sim_params`: A struct with simulation parameters.
"""
function create_uciwweihr_model_params(; 
    E_init_sd::Float64=50.0, E_init_mean::Int64=200,
    I_init_sd::Float64=20.0, I_init_mean::Int64=100,
    H_init_sd::Float64=5.0, H_init_mean::Int64=20,
    gamma_sd::Float64=0.02, log_gamma_mean::Float64=log(1/4),
    nu_sd::Float64=0.02, log_nu_mean::Float64=log(1/7),
    epsilon_sd::Float64=0.02, log_epsilon_mean::Float64=log(1/5),
    rho_gene_sd::Float64=0.02, log_rho_gene_mean::Float64=log(0.011),
    tau_sd::Float64=0.02, log_tau_mean::Float64=log(0.1),
    df_shape::Float64=2.0, df_scale::Float64=10.0,
    sigma_hosp_sd::Float64=50.0, sigma_hosp_mean::Float64=500.0,
    Rt_init_sd::Float64=0.3, Rt_init_mean::Float64=0.2,
    sigma_Rt_sd::Float64=0.2, sigma_Rt_mean::Float64=-3.0,
    w_init_sd::Float64=0.04, w_init_mean::Float64=logit(0.35),
    sigma_w_sd::Float64=0.2, sigma_w_mean::Float64=-3.5
    )

    
    return uciwweihr_model_params(
        E_init_sd, E_init_mean,
        I_init_sd, I_init_mean,
        H_init_sd, H_init_mean,
        gamma_sd, log_gamma_mean,
        nu_sd, log_nu_mean,
        epsilon_sd, log_epsilon_mean,
        rho_gene_sd, log_rho_gene_mean,
        tau_sd, log_tau_mean,
        df_shape, df_scale,
        sigma_hosp_sd, sigma_hosp_mean,
        Rt_init_sd, Rt_init_mean,
        sigma_Rt_sd, sigma_Rt_mean,
        w_init_sd, w_init_mean,
        sigma_w_sd, sigma_w_mean
    )
end