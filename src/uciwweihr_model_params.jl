"""
    uciwweihr_model_params1

Struct for holding parameters used in the UCIWWEIHR ODE compartmental model.  Use `create_uciwweihr_model_params1` to create an instance of this struct.  With hard coded sigma_wastewater and sigma_hosp.
This will use a time-varying hospitalization probability.

# Fields
- `E_init_sd::Float64=50.0`: Standard deviation for the initial number of exposed individuals.
- `log_E_init_mean::Int64=200`: Mean for the initial number of exposed individuals, on log scale.
- `I_init_sd::Float64=20.0`: Standard deviation for the initial number of infected individuals.
- `log_I_init_mean::Int64=100`: Mean for the initial number of infected individuals, on log scale.
- `H_init_sd::Float64=5.0`: Standard deviation for the initial number of hospitalized individuals.
- `log_H_init_mean::Int64=20`: Mean for the initial number of hospitalized individuals, on log scale.
- `gamma_sd::Float64=0.02`: Standard deviation for the rate of incubation.
- `log_gamma_mean::Float64=log(1/4)`: Mean for the rate of incubation on log scale.
- `nu_sd::Float64=0.02`: Standard deviation for the rate of leaving the infected compartment.
- `log_nu_mean::Float64=log(1/7)`: Mean for the rate of leaving the infected compartment on the log scale.
- `epsilon_sd::Float64=0.02`: Standard deviation for the rate of hospitalization recovery.
- `log_epsilon_mean::Float64=log(1/5)`: Mean for the rate of hospitalization recovery on the log scale.
- `rho_gene_sd::Float64=0.02`: Standard deviation for the rho prior.
- `log_rho_gene_mean::Float64=log(0.011)`: Mean for the row prior on log scale.
- `sigma_ww::Float64=log(0.1)`: Standard deviation for the normal distribution for wastewater data.  Not infered.
- `sigma_hosp::Float64=500.0`: Standard deviation for the negative binomial distribution for hospital data.  Not infered.
- `Rt_init_sd::Float64=0.3`: Standard deviation for the initial value of the time-varying reproduction number.
- `Rt_init_mean::Float64=0.2`: Mean for the initial value of the time-varying reproduction number.
- `sigma_Rt_sd::Float64=0.2`: Standard deviation for normal prior of log time-varying reproduction number standard deviation.
- `sigma_Rt_mean::Float64=-3.0`: Mean for normal prior of log time-varying reproduction number standard deviation.
- `w_init_sd::Float64=0.1`: Standard deviation for the initial value of the time-varying hospitalization rate.
- `w_init_mean::Float64=log(0.35)`: Mean for the initial value of the time-varying hospitalization rate.
- `sigma_w_sd::Float64=0.2`: Standard deviation for normal prior of log time-varying hospitalization rate standard deviation.
- `sigma_w_mean::Float64=-3.5`: Mean for normal prior of time-varying hospitalization rate standard deviation.
"""
struct uciwweihr_model_params1
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

    sigma_wastewater::Float64
    sigma_hosp::Float64

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
    uciwweihr_model_params2

Struct for holding parameters used in the UCIWWEIHR ODE compartmental model.  Use `create_uciwweihr_model_params2` to create an instance of this struct.  With prior on sigma_wastewater and sigma_hosp.
This will use a time-varying hospitalization probability.

# Fields
- `E_init_sd::Float64=50.0`: Standard deviation for the initial number of exposed individuals.
- `log_E_init_mean::Int64=200`: Mean for the initial number of exposed individuals, on log scale.
- `I_init_sd::Float64=20.0`: Standard deviation for the initial number of infected individuals.
- `log_I_init_mean::Int64=100`: Mean for the initial number of infected individuals, on log scale.
- `H_init_sd::Float64=5.0`: Standard deviation for the initial number of hospitalized individuals.
- `log_H_init_mean::Int64=20`: Mean for the initial number of hospitalized individuals, on log scale.
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
- `sigma_hosp_mean::Float64=nothing`: Mean for the normal prior of the log standard deviation of the hospital data. If `nothing`, the sigma_hosp is used.
- `Rt_init_sd::Float64=0.3`: Standard deviation for the initial value of the time-varying reproduction number.
- `Rt_init_mean::Float64=0.2`: Mean for the initial value of the time-varying reproduction number.
- `sigma_Rt_sd::Float64=0.2`: Standard deviation for normal prior of log time-varying reproduction number standard deviation.
- `sigma_Rt_mean::Float64=-3.0`: Mean for normal prior of log time-varying reproduction number standard deviation.
- `w_init_sd::Float64=0.1`: Standard deviation for the initial value of the time-varying hospitalization rate.
- `w_init_mean::Float64=log(0.35)`: Mean for the initial value of the time-varying hospitalization rate.
- `sigma_w_sd::Float64=0.2`: Standard deviation for normal prior of log time-varying hospitalization rate standard deviation.
- `sigma_w_mean::Float64=-3.5`: Mean for normal prior of time-varying hospitalization rate standard deviation.
"""
struct uciwweihr_model_params2
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
    uciwweihr_model_params3

Struct for holding parameters used in the UCIWWEIHR ODE compartmental model.  Use `create_uciwweihr_model_params1` to create an instance of this struct.  With hard coded sigma_wastewater and sigma_hosp.
This not will use a time-varying hospitalization probability.

# Fields
- `E_init_sd::Float64=50.0`: Standard deviation for the initial number of exposed individuals.
- `log_E_init_mean::Int64=200`: Mean for the initial number of exposed individuals, on log scale.
- `I_init_sd::Float64=20.0`: Standard deviation for the initial number of infected individuals.
- `log_I_init_mean::Int64=100`: Mean for the initial number of infected individuals, on log scale.
- `H_init_sd::Float64=5.0`: Standard deviation for the initial number of hospitalized individuals.
- `log_H_init_mean::Int64=20`: Mean for the initial number of hospitalized individuals, on log scale.
- `gamma_sd::Float64=0.02`: Standard deviation for the rate of incubation.
- `log_gamma_mean::Float64=log(1/4)`: Mean for the rate of incubation on log scale.
- `nu_sd::Float64=0.02`: Standard deviation for the rate of leaving the infected compartment.
- `log_nu_mean::Float64=log(1/7)`: Mean for the rate of leaving the infected compartment on the log scale.
- `epsilon_sd::Float64=0.02`: Standard deviation for the rate of hospitalization recovery.
- `log_epsilon_mean::Float64=log(1/5)`: Mean for the rate of hospitalization recovery on the log scale.
- `rho_gene_sd::Float64=0.02`: Standard deviation for the rho prior.
- `log_rho_gene_mean::Float64=log(0.011)`: Mean for the row prior on log scale.
- `sigma_ww::Float64=log(0.1)`: Standard deviation for the normal distribution for wastewater data.  Not infered.
- `sigma_hosp::Float64=500.0`: Standard deviation for the negative binomial distribution for hospital data.  Not infered.
- `Rt_init_sd::Float64=0.3`: Standard deviation for the initial value of the time-varying reproduction number.
- `Rt_init_mean::Float64=0.2`: Mean for the initial value of the time-varying reproduction number.
- `sigma_Rt_sd::Float64=0.2`: Standard deviation for normal prior of log time-varying reproduction number standard deviation.
- `sigma_Rt_mean::Float64=-3.0`: Mean for normal prior of log time-varying reproduction number standard deviation.
- `w_sd::Float64=0.1`: Standard deviation for the initial value of the time-varying hospitalization probability.
- `logit_w_mean::Float64=log(0.35)`: Mean for the initial value of the time-varying hospitalization probability. On the logit scale.
"""
struct uciwweihr_model_params3
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

    sigma_wastewater::Float64
    sigma_hosp::Float64

    Rt_init_sd::Float64
    Rt_init_mean::Float64
    sigma_Rt_sd::Float64
    sigma_Rt_mean::Float64
    w_sd::Float64
    logit_w_mean::Float64
end

"""
    uciwweihr_model_params4

Struct for holding parameters used in the UCIWWEIHR ODE compartmental model.  Use `create_uciwweihr_model_params2` to create an instance of this struct.  With prior on sigma_wastewater and sigma_hosp.
This will not use a time-varying hospitalization probability.

# Fields
- `E_init_sd::Float64=50.0`: Standard deviation for the initial number of exposed individuals.
- `log_E_init_mean::Int64=200`: Mean for the initial number of exposed individuals, on log scale.
- `I_init_sd::Float64=20.0`: Standard deviation for the initial number of infected individuals.
- `log_I_init_mean::Int64=100`: Mean for the initial number of infected individuals, on log scale.
- `H_init_sd::Float64=5.0`: Standard deviation for the initial number of hospitalized individuals.
- `log_H_init_mean::Int64=20`: Mean for the initial number of hospitalized individuals, on log scale.
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
- `sigma_hosp_mean::Float64=nothing`: Mean for the normal prior of the log standard deviation of the hospital data. If `nothing`, the sigma_hosp is used.
- `Rt_init_sd::Float64=0.3`: Standard deviation for the initial value of the time-varying reproduction number.
- `Rt_init_mean::Float64=0.2`: Mean for the initial value of the time-varying reproduction number.
- `sigma_Rt_sd::Float64=0.2`: Standard deviation for normal prior of log time-varying reproduction number standard deviation.
- `sigma_Rt_mean::Float64=-3.0`: Mean for normal prior of log time-varying reproduction number standard deviation.
- `w_sd::Float64=0.1`: Standard deviation for the initial value of the time-varying hospitalization probability.
- `logit_w_mean::Float64=log(0.35)`: Mean for the initial value of the time-varying hospitalization probability. On the logit scale.
"""
struct uciwweihr_model_params4
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
    sigma_hosp_mean::Float64

    Rt_init_sd::Float64
    Rt_init_mean::Float64
    sigma_Rt_sd::Float64
    sigma_Rt_mean::Float64
    w_sd::Float64
    logit_w_mean::Float64
end






"""
    create_uciwweihr_model_params1(; kwargs...) or create_uciwweihr_model_params2(; kwargs...) or create_uciwweihr_model_params3(; kwargs...) or create_uciwweihr_model_params4(; kwargs...)

Creates a `uciwweihr_sim_params1` or `uciwweihr_sim_params2` struct with the option to either hard code sigma_wastewater and sigma_hosp or provide a prior with time-varying hospitalization probability.  
Can also be used to create a `uciwweihr_sim_params3` or `uciwweihr_sim_params4` struct with the option to either hard code sigma_wastewater and sigma_hosp or provide a prior without time-varying hospitalization probability.

# Arguments
- `kwargs...`: Named arguments corresponding to the fields in `uciwweihr_sim_params1` or `uciwweihr_sim_params2`.

# Returns
- `params::uciwweihr_sim_params1` or `params::uciwweihr_sim_params2`: A struct with simulation parameters.
"""
function create_uciwweihr_model_params2(; 
    E_init_sd::Float64=0.2, log_E_init_mean::Float64=log(200),
    I_init_sd::Float64=0.2, log_I_init_mean::Float64=log(100),
    H_init_sd::Float64=0.2, log_H_init_mean::Float64=log(20),
    gamma_sd::Float64=0.02, log_gamma_mean::Float64=log(1/4),
    nu_sd::Float64=0.02, log_nu_mean::Float64=log(1/7),
    epsilon_sd::Float64=0.02, log_epsilon_mean::Float64=log(1/5),
    rho_gene_sd::Float64=0.02, log_rho_gene_mean::Float64=log(0.011),
    
    sigma_ww_sd::Float64=0.02, log_sigma_ww_mean::Float64=log(0.1),
    sigma_hosp_sd::Float64=50.0, sigma_hosp_mean::Float64=500.0,

    Rt_init_sd::Float64=0.3, Rt_init_mean::Float64=0.2,
    sigma_Rt_sd::Float64=0.2, sigma_Rt_mean::Float64=-3.0,
    w_init_sd::Float64=0.04, w_init_mean::Float64=logit(0.35),
    sigma_w_sd::Float64=0.2, sigma_w_mean::Float64=-3.5,
    message::Bool=true
    )
    if message
        println("Using prior for sigma_ww and sigma_hosp / with time-varying hospitalization probability!!!") 
    end
    
    return uciwweihr_model_params2(
        E_init_sd, log_E_init_mean,
        I_init_sd, log_I_init_mean,
        H_init_sd, log_H_init_mean,
        gamma_sd, log_gamma_mean,
        nu_sd, log_nu_mean,
        epsilon_sd, log_epsilon_mean,
        rho_gene_sd, log_rho_gene_mean,

        sigma_ww_sd, log_sigma_ww_mean,
        sigma_hosp_sd, sigma_hosp_mean,

        Rt_init_sd, Rt_init_mean,
        sigma_Rt_sd, sigma_Rt_mean,
        w_init_sd, w_init_mean,
        sigma_w_sd, sigma_w_mean
    )
end


function create_uciwweihr_model_params1(; 
    E_init_sd::Float64=0.2, log_E_init_mean::Float64=log(200),
    I_init_sd::Float64=0.2, log_I_init_mean::Float64=log(100),
    H_init_sd::Float64=0.2, log_H_init_mean::Float64=log(20),
    gamma_sd::Float64=0.02, log_gamma_mean::Float64=log(1/4),
    nu_sd::Float64=0.02, log_nu_mean::Float64=log(1/7),
    epsilon_sd::Float64=0.02, log_epsilon_mean::Float64=log(1/5),
    rho_gene_sd::Float64=0.02, log_rho_gene_mean::Float64=log(0.011),

    sigma_wastewater::Float64=0.1,
    sigma_hosp::Float64=500.0,

    Rt_init_sd::Float64=0.3, Rt_init_mean::Float64=0.2,
    sigma_Rt_sd::Float64=0.2, sigma_Rt_mean::Float64=-3.0,
    w_init_sd::Float64=0.04, w_init_mean::Float64=logit(0.35),
    sigma_w_sd::Float64=0.2, sigma_w_mean::Float64=-3.5,
    message::Bool=true
    )
    if message
        println("Using hard coded for sigma_ww and sigma_hosp / with time-varying hospitalization probability!!!") 
    end
    
    return uciwweihr_model_params1(
        E_init_sd, log_E_init_mean,
        I_init_sd, log_I_init_mean,
        H_init_sd, log_H_init_mean,
        gamma_sd, log_gamma_mean,
        nu_sd, log_nu_mean,
        epsilon_sd, log_epsilon_mean,
        rho_gene_sd, log_rho_gene_mean,

        sigma_wastewater, sigma_hosp,

        Rt_init_sd, Rt_init_mean,
        sigma_Rt_sd, sigma_Rt_mean,
        w_init_sd, w_init_mean,
        sigma_w_sd, sigma_w_mean
    )
end



function create_uciwweihr_model_params4(; 
    E_init_sd::Float64=0.2, log_E_init_mean::Float64=log(200),
    I_init_sd::Float64=0.2, log_I_init_mean::Float64=log(100),
    H_init_sd::Float64=0.2, log_H_init_mean::Float64=log(20),
    gamma_sd::Float64=0.02, log_gamma_mean::Float64=log(1/4),
    nu_sd::Float64=0.02, log_nu_mean::Float64=log(1/7),
    epsilon_sd::Float64=0.02, log_epsilon_mean::Float64=log(1/5),
    rho_gene_sd::Float64=0.02, log_rho_gene_mean::Float64=log(0.011),
    
    sigma_ww_sd::Float64=0.02, log_sigma_ww_mean::Float64=log(0.1),
    sigma_hosp_sd::Float64=50.0, sigma_hosp_mean::Float64=500.0,

    Rt_init_sd::Float64=0.3, Rt_init_mean::Float64=0.2,
    sigma_Rt_sd::Float64=0.2, sigma_Rt_mean::Float64=-3.0,
    w_sd::Float64=0.04, logit_w_mean::Float64=logit(0.35),
    message::Bool=true
    )
    if message
        println("Using prior for sigma_ww and sigma_hosp / without time-varying hospitalization probability!!!") 
    end
    
    return uciwweihr_model_params4(
        E_init_sd, log_E_init_mean,
        I_init_sd, log_I_init_mean,
        H_init_sd, log_H_init_mean,
        gamma_sd, log_gamma_mean,
        nu_sd, log_nu_mean,
        epsilon_sd, log_epsilon_mean,
        rho_gene_sd, log_rho_gene_mean,

        sigma_ww_sd, log_sigma_ww_mean,
        sigma_hosp_sd, sigma_hosp_mean,

        Rt_init_sd, Rt_init_mean,
        sigma_Rt_sd, sigma_Rt_mean,
        w_sd, logit_w_mean
    )
end


function create_uciwweihr_model_params3(; 
    E_init_sd::Float64=0.2, log_E_init_mean::Float64=log(200),
    I_init_sd::Float64=0.2, log_I_init_mean::Float64=log(100),
    H_init_sd::Float64=0.2, log_H_init_mean::Float64=log(20),
    gamma_sd::Float64=0.02, log_gamma_mean::Float64=log(1/4),
    nu_sd::Float64=0.02, log_nu_mean::Float64=log(1/7),
    epsilon_sd::Float64=0.02, log_epsilon_mean::Float64=log(1/5),
    rho_gene_sd::Float64=0.02, log_rho_gene_mean::Float64=log(0.011),

    sigma_wastewater::Float64=0.1,
    sigma_hosp::Float64=500.0,

    Rt_init_sd::Float64=0.3, Rt_init_mean::Float64=0.2,
    sigma_Rt_sd::Float64=0.2, sigma_Rt_mean::Float64=-3.0,
    w_sd::Float64=0.04, logit_w_mean::Float64=logit(0.35),
    message::Bool=true
    )
    if message
        println("Using hard coded for sigma_ww and sigma_hosp / without time-varying hospitalization probability!!!") 
    end
    
    return uciwweihr_model_params3(
        E_init_sd, log_E_init_mean,
        I_init_sd, log_I_init_mean,
        H_init_sd, log_H_init_mean,
        gamma_sd, log_gamma_mean,
        nu_sd, log_nu_mean,
        epsilon_sd, log_epsilon_mean,
        rho_gene_sd, log_rho_gene_mean,

        sigma_wastewater, sigma_hosp,

        Rt_init_sd, Rt_init_mean,
        sigma_Rt_sd, sigma_Rt_mean,
        w_sd, logit_w_mean
    )
end