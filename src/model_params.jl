"""
    model_params_time_var_hosp

Struct for holding parameters used in the UCIWWEIHR ODE compartmental model.  Use `create_model_params_time_var_hosp` to create an instance of this struct.
This will use a time-varying hospitalization probability.

# Fields
- `E_init_sd::Float64=50.0`: Standard deviation for the initial number of exposed individuals.
- `log_E_init_mean::Int64=200`: Mean for the initial number of exposed individuals, on log scale.
- `I_init_sd::Float64=20.0`: Standard deviation for the initial number of infected individuals.
- `log_I_init_mean::Int64=100`: Mean for the initial number of infected individuals, on log scale.
- `H_init_sd::Float64=5.0`: Standard deviation for the initial number of hospitalized individuals.
- `log_H_init_mean::Int64=20`: Mean for the initial number of hospitalized individuals, on log scale.
- `CH_init::Float64=5.0`: Initial number of cumulative hospitalized individuals.
- `log_CH_init_mean::Int64=20`: Mean for the initial number of cumulative hospitalized individuals, on log scale.
- `gamma_sd::Float64=0.02`: Standard deviation for the rate of incubation.
- `log_gamma_mean::Float64=log(1/4)`: Mean for the rate of incubation on log scale.
- `nu_sd::Float64=0.02`: Standard deviation for the rate of leaving the infected compartment.
- `log_nu_mean::Float64=log(1/7)`: Mean for the rate of leaving the infected compartment on the log scale.
- `epsilon_sd::Float64=0.02`: Standard deviation for the rate of hospitalization recovery.
- `log_epsilon_mean::Float64=log(1/5)`: Mean for the rate of hospitalization recovery on the log scale.
- `rho_gene_sd::Float64=0.02`: Standard deviation for the rho prior.
- `log_rho_gene_mean::Float64=log(0.011)`: Mean for the row prior on log scale.
- `sigma_ww_sd::Float64=nothing`: Standard deviation for the normal prior of the log standard deviation of the wastewater data. If `nothing`, the sigma_ww is used.
- `log_sigma_ww_mean::Float64=nothing`: Mean for the normal prior of the log standard deviation of the wastewater data. If `nothing`, the sigma_ww is used.
- `sigma_hosp_sd::Float64=nothing`: Standard deviation for the normal prior of the log standard deviation of the hospital data. If `nothing`, the sigma_hosp is used.
- `log_sigma_hosp_mean::Float64=nothing`: Mean for the normal prior of the log standard deviation of the hospital data. If `nothing`, the sigma_hosp is used.
- `Rt_init_sd::Float64=0.3`: Standard deviation for the initial value of the time-varying reproduction number.
- `Rt_init_mean::Float64=0.2`: Mean for the initial value of the time-varying reproduction number.
- `sigma_Rt_sd::Float64=0.2`: Standard deviation for normal prior of log time-varying reproduction number standard deviation.
- `sigma_Rt_mean::Float64=-3.0`: Mean for normal prior of log time-varying reproduction number standard deviation.
- `w_init_sd::Float64=0.1`: Standard deviation for the initial value of the time-varying hospitalization rate.
- `w_init_mean::Float64=log(0.35)`: Mean for the initial value of the time-varying hospitalization rate.
- `sigma_w_sd::Float64=0.2`: Standard deviation for normal prior of log time-varying hospitalization rate standard deviation.
- `sigma_w_mean::Float64=-3.5`: Mean for normal prior of time-varying hospitalization rate standard deviation.
"""
struct model_params_time_var_hosp_prev
    E_init_sd::Float64
    log_E_init_mean::Float64
    I_init_sd::Float64
    log_I_init_mean::Float64
    H_init_sd::Float64
    log_H_init_mean::Float64
    gamma_sd::Float64
    log_gamma_mean::Float64
    nu_sd::Float64
    log_nu_mean::Float64
    epsilon_sd::Float64
    log_epsilon_mean::Float64
    rho_gene_sd::Float64
    log_rho_gene_mean::Float64

    sigma_ww_sd::Float64
    log_sigma_ww_mean::Float64
    sigma_hosp_sd::Float64
    log_sigma_hosp_mean::Float64

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
    create_model_params_time_var_hosp(; kwargs...) or create_model_params_non_time_var_hosp(; kwargs...)

Creates a `model_params_time_var_hosp` or `model_params2` struct with the option to either have time-varying hospitalization probability or not.
# Arguments
- `kwargs...`: Named arguments corresponding to the fields in `model_params1` or `model_params2`.

# Returns
- `params::model_params_time_var_hosp` or `params::model_params_non_time_var_hosp`: A struct with simulation parameters.
"""
function create_model_params_time_var_hosp_prev(
    E_init_sd::Float64, log_E_init_mean::Float64,
    I_init_sd::Float64, log_I_init_mean::Float64,
    H_init_sd::Float64, log_H_init_mean::Float64,
    gamma_sd::Float64, log_gamma_mean::Float64,
    nu_sd::Float64, log_nu_mean::Float64,
    epsilon_sd::Float64, log_epsilon_mean::Float64,
    rho_gene_sd::Float64, log_rho_gene_mean::Float64,

    sigma_ww_sd::Float64, log_sigma_ww_mean::Float64,
    sigma_hosp_sd::Float64, log_sigma_hosp_mean::Float64,

    Rt_init_sd::Float64, Rt_init_mean::Float64,
    sigma_Rt_sd::Float64, sigma_Rt_mean::Float64,
    w_init_sd::Float64, w_init_mean::Float64,
    sigma_w_sd::Float64, sigma_w_mean::Float64,
    message::Bool;
    )
    if message
        println("Using time-varying hospitalization probability - Prevalence Model parameters!!!") 
    end

    return model_params_time_var_hosp_prev(
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
        sigma_w_sd, sigma_w_mean
    )
end


struct model_params_time_var_hosp_inc
    E_init_sd::Float64
    log_E_init_mean::Float64
    I_init_sd::Float64
    log_I_init_mean::Float64
    gamma_sd::Float64
    log_gamma_mean::Float64
    nu_sd::Float64
    log_nu_mean::Float64
    rho_gene_sd::Float64
    log_rho_gene_mean::Float64

    sigma_ww_sd::Float64
    log_sigma_ww_mean::Float64
    sigma_hosp_sd::Float64
    log_sigma_hosp_mean::Float64

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
    create_model_params_time_var_hosp(; kwargs...) or create_model_params_non_time_var_hosp(; kwargs...)

Creates a `model_params_time_var_hosp` or `model_params2` struct with the option to either have time-varying hospitalization probability or not.
# Arguments
- `kwargs...`: Named arguments corresponding to the fields in `model_params1` or `model_params2`.

# Returns
- `params::model_params_time_var_hosp` or `params::model_params_non_time_var_hosp`: A struct with simulation parameters.
"""
function create_model_params_time_var_hosp_inc(
    E_init_sd::Float64, log_E_init_mean::Float64,
    I_init_sd::Float64, log_I_init_mean::Float64,
    gamma_sd::Float64, log_gamma_mean::Float64,
    nu_sd::Float64, log_nu_mean::Float64,
    rho_gene_sd::Float64, log_rho_gene_mean::Float64,

    sigma_ww_sd::Float64, log_sigma_ww_mean::Float64,
    sigma_hosp_sd::Float64, log_sigma_hosp_mean::Float64,

    Rt_init_sd::Float64, Rt_init_mean::Float64,
    sigma_Rt_sd::Float64, sigma_Rt_mean::Float64,
    w_init_sd::Float64, w_init_mean::Float64,
    sigma_w_sd::Float64, sigma_w_mean::Float64,
    message::Bool;
    )
    if message
        println("Using time-varying hospitalization probability - Incidence Model parameters!!!") 
    end

    return model_params_time_var_hosp_inc(
        E_init_sd, log_E_init_mean,
        I_init_sd, log_I_init_mean,
        gamma_sd, log_gamma_mean,
        nu_sd, log_nu_mean,
        rho_gene_sd, log_rho_gene_mean,

        sigma_ww_sd, log_sigma_ww_mean,
        sigma_hosp_sd, log_sigma_hosp_mean,

        Rt_init_sd, Rt_init_mean,
        sigma_Rt_sd, sigma_Rt_mean,
        w_init_sd, w_init_mean,
        sigma_w_sd, sigma_w_mean
    )
end


"""
    model_params_non_time_var_hosp_no_ww

Struct for holding parameters used in the UCIWWEIHR ODE compartmental model without wastewater.  Use `create_model_params_non_time_var_hosp_no_ww` to create an instance of this struct.
This will not use a time-varying hospitalization probability.

# Fields
- `E_init_sd::Float64=50.0`: Standard deviation for the initial number of exposed individuals.
- `log_E_init_mean::Int64=200`: Mean for the initial number of exposed individuals, on log scale.
- `I_init_sd::Float64=20.0`: Standard deviation for the initial number of infected individuals.
- `log_I_init_mean::Int64=100`: Mean for the initial number of infected individuals, on log scale.
- `H_init_sd::Float64=5.0`: Standard deviation for the initial number of hospitalized individuals.
- `log_H_init_mean::Int64=20`: Mean for the initial number of hospitalized individuals, on log scale.
- `CH_init::Float64=5.0`: Initial number of cumulative hospitalized individuals.
- `gamma_sd::Float64=0.02`: Standard deviation for the rate of incubation.
- `log_gamma_mean::Float64=log(1/4)`: Mean for the rate of incubation on log scale.
- `nu_sd::Float64=0.02`: Standard deviation for the rate of leaving the infected compartment.
- `log_nu_mean::Float64=log(1/7)`: Mean for the rate of leaving the infected compartment on the log scale.
- `epsilon_sd::Float64=0.02`: Standard deviation for the rate of hospitalization recovery.
- `log_epsilon_mean::Float64=log(1/5)`: Mean for the rate of hospitalization recovery on the log scale.
- `sigma_hosp_sd::Float64=nothing`: Standard deviation for the normal prior of the log standard deviation of the hospital data. If `nothing`, the sigma_hosp is used.
- `log_sigma_hosp_mean::Float64=nothing`: Mean for the normal prior of the log standard deviation of the hospital data. If `nothing`, the sigma_hosp is used.
- `Rt_init_sd::Float64=0.3`: Standard deviation for the initial value of the time-varying reproduction number.
- `Rt_init_mean::Float64=0.2`: Mean for the initial value of the time-varying reproduction number.
- `sigma_Rt_sd::Float64=0.2`: Standard deviation for normal prior of log time-varying reproduction number standard deviation.
- `sigma_Rt_mean::Float64=-3.0`: Mean for normal prior of log time-varying reproduction number standard deviation.
- `w_sd::Float64=0.1`: Standard deviation for the initial value of the time-varying hospitalization probability.
- `logit_w_mean::Float64=log(0.35)`: Mean for the initial value of the time-varying hospitalization probability. On the logit scale.
"""
struct model_params_non_time_var_hosp_no_ww
    E_init_sd::Float64
    log_E_init_mean::Float64
    I_init_sd::Float64
    log_I_init_mean::Float64
    H_init_sd::Float64
    log_H_init_mean::Float64
    CH_init::Float64
    gamma_sd::Float64
    log_gamma_mean::Float64
    nu_sd::Float64
    log_nu_mean::Float64
    epsilon_sd::Float64
    log_epsilon_mean::Float64

    sigma_hosp_sd::Float64
    log_sigma_hosp_mean::Float64

    Rt_init_sd::Float64
    Rt_init_mean::Float64
    sigma_Rt_sd::Float64
    sigma_Rt_mean::Float64
    w_sd::Float64
    logit_w_mean::Float64
end


function create_model_params_non_time_var_hosp(
    E_init_sd::Float64, log_E_init_mean::Float64,
    I_init_sd::Float64, log_I_init_mean::Float64,
    H_init_sd::Float64, log_H_init_mean::Float64,
    CH_init::Float64,
    gamma_sd::Float64, log_gamma_mean::Float64,
    nu_sd::Float64, log_nu_mean::Float64,
    epsilon_sd::Float64, log_epsilon_mean::Float64,
    
    sigma_hosp_sd::Float64, log_sigma_hosp_mean::Float64,

    Rt_init_sd::Float64, Rt_init_mean::Float64,
    sigma_Rt_sd::Float64, sigma_Rt_mean::Float64,
    w_sd::Float64, logit_w_mean::Float64,
    message::Bool;
    )
    if message
        println("Using model without time-varying hospitalization probability w/out wastewater parameters!!!") 
    end

    return model_params_non_time_var_hosp_no_ww(
        E_init_sd, log_E_init_mean,
        I_init_sd, log_I_init_mean,
        H_init_sd, log_H_init_mean,
        CH_init,
        gamma_sd, log_gamma_mean,
        nu_sd, log_nu_mean,
        epsilon_sd, log_epsilon_mean,

        sigma_hosp_sd, log_sigma_hosp_mean,

        Rt_init_sd, Rt_init_mean,
        sigma_Rt_sd, sigma_Rt_mean,
        w_sd, logit_w_mean
    )
end
