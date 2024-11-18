"""
    uciwweihr_model_params

Struct for holding parameters used in the UCIWWEIHR ODE compartmental model.  Use `create_uciwweihr_model_params` to create an instance of this struct.

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
- `sigma_ww::Union{Float64, Nothing}=log(0.1)`: Standard deviation for the normal distribution for wastewater data.  Not infered.
- `sigma_hosp::Union{Float64, Nothing}=500.0`: Standard deviation for the negative binomial distribution for hospital data.  Not infered.
- `sigma_ww_sd::Union{Float64, Nothing}=nothing`: Standard deviation for the normal prior of the log standard deviation of the wastewater data. If `nothing`, the sigma_ww is used.
- `log_sigma_ww_mean::Union{Float64, Nothing}=nothing`: Mean for the normal prior of the log standard deviation of the wastewater data. If `nothing`, the sigma_ww is used.
- `sigma_hosp_sd::Union{Float64, Nothing}=nothing`: Standard deviation for the normal prior of the log standard deviation of the hospital data. If `nothing`, the sigma_hosp is used.
- `sigma_hosp_mean::Union{Float64, Nothing}=nothing`: Mean for the normal prior of the log standard deviation of the hospital data. If `nothing`, the sigma_hosp is used.
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

    sigma_ww::Union{Float64, Nothing}
    sigma_hosp::Union{Float64, Nothing}

    sigma_ww_sd::Union{Float64, Nothing}
    log_sigma_ww_mean::Union{Float64, Nothing}
    sigma_hosp_sd::Union{Float64, Nothing}
    sigma_hosp_mean::Union{Float64, Nothing}

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
    E_init_sd::Float64=50.0, log_E_init_mean::Float64=log(200),
    I_init_sd::Float64=20.0, log_I_init_mean::Float64=log(100),
    H_init_sd::Float64=5.0, log_H_init_mean::Float64=log(20),
    gamma_sd::Float64=0.02, log_gamma_mean::Float64=log(1/4),
    nu_sd::Float64=0.02, log_nu_mean::Float64=log(1/7),
    epsilon_sd::Float64=0.02, log_epsilon_mean::Float64=log(1/5),
    rho_gene_sd::Float64=0.02, log_rho_gene_mean::Float64=log(0.011),

    sigma_wastewater::Union{Float64, Nothing}=0.1,
    sigma_hosp::Union{Float64, Nothing}=500.0,
    
    #sigma_ww_sd::Float64=0.02, log_sigma_ww_mean::Float64=log(0.1),
    #sigma_hosp_sd::Float64=50.0, sigma_hosp_mean::Float64=500.0,
    sigma_ww_sd::Union{Float64, Nothing}=nothing, log_sigma_ww_mean::Union{Float64, Nothing}=nothing,
    sigma_hosp_sd::Union{Float64, Nothing}=nothing, sigma_hosp_mean::Union{Float64, Nothing}=nothing,


    Rt_init_sd::Float64=0.3, Rt_init_mean::Float64=0.2,
    sigma_Rt_sd::Float64=0.2, sigma_Rt_mean::Float64=-3.0,
    w_init_sd::Float64=0.04, w_init_mean::Float64=logit(0.35),
    sigma_w_sd::Float64=0.2, sigma_w_mean::Float64=-3.5
    )

    if sigma_ww_sd == nothing
        println("sigma_ww will be fixed")
    end
    if sigma_hosp_sd == nothing
        println("sigma_hosp will be fixed")
    end
    
    return uciwweihr_model_params(
        E_init_sd, log_E_init_mean,
        I_init_sd, log_I_init_mean,
        H_init_sd, log_H_init_mean,
        gamma_sd, log_gamma_mean,
        nu_sd, log_nu_mean,
        epsilon_sd, log_epsilon_mean,
        rho_gene_sd, log_rho_gene_mean,
        sigma_wastewater, sigma_hosp,

        sigma_ww_sd, log_sigma_ww_mean,
        sigma_hosp_sd, sigma_hosp_mean,

        Rt_init_sd, Rt_init_mean,
        sigma_Rt_sd, sigma_Rt_mean,
        w_init_sd, w_init_mean,
        sigma_w_sd, sigma_w_mean
    )
end
