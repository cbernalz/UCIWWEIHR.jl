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
include("generate_simulation_data_agent.jl")
include("eihr_ode.jl")
include("negativebinomial2.jl")
include("generalizedtdist.jl")
include("uciwweihr_model_params.jl")
include("uciwweihr_init_param.jl")
include("uciwweihr_likelihood_helpers.jl")
include("uciwweihr_model.jl")
include("uciwweihr_fit.jl")
include("uciwweihr_gq_pp.jl")
include("optimize_many_MAP.jl")
include("helper_functions.jl")
include("mcmcdiags_vis.jl")
include("time_varying_param_vis.jl")
include("non_time_varying_param_vis.jl")
include("predictive_param_vis.jl")
include("ode_solution_vis.jl")
include("uciwweihr_visualizer.jl")

export eihr_ode
export uciwweihr_sim_params
export create_uciwweihr_sim_params
export generate_random_walk
export generate_logit_normal_random_walk
export generate_simulation_data_uciwweihr
export generate_simulation_data_agent
export NegativeBinomial2
export GeneralizedTDist
export uciwweihr_model_params1
export uciwweihr_model_params2
export uciwweihr_init_param
export uciwweihr_likelihood_helpers
export create_uciwweihr_model_params1
export create_uciwweihr_model_params2
export uciwweihr_model
export uciwweihr_fit
export uciwweihr_gq_pp
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