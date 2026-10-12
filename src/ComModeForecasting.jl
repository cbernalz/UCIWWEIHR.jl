module ComModeForecasting   

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
#using ForwardDiff
using Logging
using CSV
using DataFrames
#using DifferentialEquations
using StatsBase
using Plots
using Printf
using Colors

include("generate_simulation_data_ComModeForecasting.jl")
include("eihr_ode.jl")
include("negativebinomial2.jl")
include("generalizedtdist.jl")
include("model_params.jl")
include("likelihood_helpers.jl")
include("ComModeForecasting_model.jl")
include("fit.jl")
include("generate_pq_pp.jl")
include("optimize_many_MAP.jl")
include("helper_functions.jl")


export eihr_ode
export eihr_ode_const_w
export eihr_ode_inc_const_w
export ComModeForecasting_sim_params
export create_ComModeForecasting_sim_params
export generate_random_walk
export generate_logit_normal_random_walk
export generate_simulation_data_ComModeForecasting
export NegativeBinomial2
export GeneralizedTDist
export model_params_time_var_hosp_prev
export model_params_time_var_hosp_inc
export model_params_non_time_var_hosp
export create_model_params_time_var_hosp_prev
export create_model_params_time_var_hosp_inc
export create_model_params_hosp_inc_prev
export create_model_params_time_var_hosp_inc_no_ww
export create_model_params_non_time_var_hosp
export likelihood_helpers
export ComModeForecasting_model
export fit
export generate_pq_pp
export optimize_many_MAP
export optimize_many_MAP2
export optimize_many_MAP2_wrapper
export ChainsCustomIndex
export save_plots_to_docs
export startswith_any
export generate_colors
export calculate_quantiles
export repeat_last_n_elements
export is_time_varying_above_n
export subset_desired_ode_from_gq
export calculate_quantiles_without_chain

end