module UCIWWEIHR   

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

include("generate_simulation_data_uciwweihr.jl")
include("eihr_ode.jl")
include("negativebinomial2.jl")
include("generalizedtdist.jl")
include("model_params.jl")
include("likelihood_helpers.jl")
include("uciwweihr_model.jl")
include("fit.jl")
include("generate_pq_pp.jl")
include("optimize_many_MAP.jl")
include("helper_functions.jl")
include("mcmcdiags_vis.jl")
include("time_varying_param_vis.jl")
include("non_time_varying_param_vis.jl")
include("predictive_param_vis.jl")
include("ode_solution_vis.jl")
include("uciwweihr_visualizer.jl")

export eihr_ode
export eihr_ode_const_w
export uciwweihr_sim_params
export create_uciwweihr_sim_params
export generate_random_walk
export generate_logit_normal_random_walk
export generate_simulation_data_uciwweihr
export NegativeBinomial2
export GeneralizedTDist
export model_params_time_var_hosp_prev
export model_params_time_var_hosp_inc
export model_params_non_time_var_hosp
export create_model_params_time_var_hosp_prev
export create_model_params_time_var_hosp_inc
export create_model_params_non_time_var_hosp
export likelihood_helpers
export uciwweihr_model
export fit
export generate_pq_pp
export optimize_many_MAP
export optimize_many_MAP2
export optimize_many_MAP2_wrapper
export uciwweihr_visualizer
export ChainsCustomIndexs
export save_plots_to_docs
export startswith_any
export generate_colors
export calculate_quantiles
export repeat_last_n_elements
export mcmcdiags_vis
export time_varying_param_vis
export non_time_varying_param_vis
export ode_solution_vis
export predictive_param_vis
export is_time_varying_above_n
export subset_desired_ode_from_gq
export calculate_quantiles_without_chain

end