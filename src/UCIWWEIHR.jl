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
using ForwardDiff
using Logging
using CSV
using DataFrames
using DifferentialEquations
using StatsBase
using Plots
using Printf
using Colors

include("generate_simulation_data_uciwweihr.jl")
include("generate_simulation_data_agent.jl")
include("eihr_ode.jl")
include("negativebinomial2.jl")
include("generalizedtdist.jl")
include("uciwweihr_model.jl")
include("uciwweihr_fit.jl")
include("uciwweihr_gq_pp.jl")
include("uciwweihr_visualizer.jl")
include("helper_functions.jl")

export eihr_ode
export generate_simulation_data_uciwweihr
export generate_simulation_data_agent
export NegativeBinomial2
export GeneralizedTDist
export uciwweihr_model
export uciwweihr_fit
export uciwweihr_gq_pp
export uciwweihr_visualizer
export ChainsCustomIndexs
export save_plots_to_docs
export startswith_any
export calculate_quantiles

end