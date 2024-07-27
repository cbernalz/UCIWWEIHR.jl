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

include("generate_simulation_data_uciwweihr.jl")
include("generate_simulation_data_agent.jl")
include("bayes_eihr_model.jl")
include("eihr_ode.jl")
include("negativebinomial2.jl")
include("generalizedtdist.jl")

export eihr_ode
export generate_simulation_data_uciwweihr
export generate_simulation_data_agent
export bayes_eihr_model
export NegativeBinomial2
export GeneralizedTDist

end