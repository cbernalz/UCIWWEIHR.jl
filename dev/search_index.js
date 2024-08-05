var documenterSearchIndex = {"docs":
[{"location":"tutorial_index/#tutorial_home","page":"TUTORIAL CONTENTS","title":"Tutorial Contents","text":"","category":"section"},{"location":"tutorial_index/","page":"TUTORIAL CONTENTS","title":"TUTORIAL CONTENTS","text":"Future Description.","category":"page"},{"location":"tutorial_index/#Contents","page":"TUTORIAL CONTENTS","title":"Contents","text":"","category":"section"},{"location":"tutorial_index/","page":"TUTORIAL CONTENTS","title":"TUTORIAL CONTENTS","text":"Getting Started\nGenerating simulated data with UCIWWEIHR ODE compartmental based model.\nGenerating simulated data with an agent based model.","category":"page"},{"location":"license/#license","page":"LICENSE","title":"UCIWWEIHR.jl Package License","text":"","category":"section"},{"location":"license/","page":"LICENSE","title":"LICENSE","text":"The UCIWWEIHR.jl package is licensed under the MIT License.","category":"page"},{"location":"license/#MIT-License","page":"LICENSE","title":"MIT License","text":"","category":"section"},{"location":"license/","page":"LICENSE","title":"LICENSE","text":"MIT License","category":"page"},{"location":"license/","page":"LICENSE","title":"LICENSE","text":"Copyright (c) 2024 Christian O. Bernal Zelaya","category":"page"},{"location":"license/","page":"LICENSE","title":"LICENSE","text":"Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"Software\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:","category":"page"},{"location":"license/","page":"LICENSE","title":"LICENSE","text":"The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.","category":"page"},{"location":"license/","page":"LICENSE","title":"LICENSE","text":"THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.","category":"page"},{"location":"tutorials/agent_based_simulation_data/#agent_based_simulation_data","page":"AGENT-BASED SIMULATION DATA","title":"Generating simulated data with an agent based model.","text":"","category":"section"},{"location":"tutorials/agent_based_simulation_data/","page":"AGENT-BASED SIMULATION DATA","title":"AGENT-BASED SIMULATION DATA","text":"This package provides a way to also simulate data using the agent based model in the future paper.  The function called generate_simulation_data_agent.jl can be used to generate synthetic data for a given population size and features.  Here we provide a demonstration using the default settings of generate_simulation_data_agent.jl :","category":"page"},{"location":"tutorials/agent_based_simulation_data/#1.-Functionality.","page":"AGENT-BASED SIMULATION DATA","title":"1. Functionality.","text":"","category":"section"},{"location":"tutorials/agent_based_simulation_data/","page":"AGENT-BASED SIMULATION DATA","title":"AGENT-BASED SIMULATION DATA","text":"using UCIWWEIHR\nusing Plots\n# Running simulation function with defaults\ndf = generate_simulation_data_agent()\nfirst(df, 5)","category":"page"},{"location":"tutorials/agent_based_simulation_data/#2.-Visualizing-SEIHR-compartments.","page":"AGENT-BASED SIMULATION DATA","title":"2. Visualizing SEIHR compartments.","text":"","category":"section"},{"location":"tutorials/agent_based_simulation_data/","page":"AGENT-BASED SIMULATION DATA","title":"AGENT-BASED SIMULATION DATA","text":"We can also use the Plots package to visualize the data generated.","category":"page"},{"location":"tutorials/agent_based_simulation_data/","page":"AGENT-BASED SIMULATION DATA","title":"AGENT-BASED SIMULATION DATA","text":"plot(df.Time, df.S, label = \"Suseptible\", \n    xlabel = \"Time\", \n    ylabel = \"Number of Individuals\", \n    title = \"Agent Based Model Simulation Results\")\nplot!(df.Time, df.E, label = \"Exposed\")\nplot!(df.Time, df.I, label = \"Infected\")\nplot!(df.Time, df.H, label = \"Hospitalized\")\nplot!(df.Time, df.R, label = \"Recovered\")","category":"page"},{"location":"tutorials/agent_based_simulation_data/#[Tutorial-Contents](@ref-tutorial_home)","page":"AGENT-BASED SIMULATION DATA","title":"Tutorial Contents","text":"","category":"section"},{"location":"tutorials/uciwweihr_model_fitting/","page":"UCIWWEIHR FITTING MODEL","title":"UCIWWEIHR FITTING MODEL","text":"using Plots, StatsPlots; gr()\nPlots.reset_defaults()\n","category":"page"},{"location":"tutorials/uciwweihr_model_fitting/#uciwwiehr_model_fitting","page":"UCIWWEIHR FITTING MODEL","title":"Generating Posterior Distribution Samples with UCIWWEIHR ODE compartmental based model.","text":"","category":"section"},{"location":"tutorials/uciwweihr_model_fitting/","page":"UCIWWEIHR FITTING MODEL","title":"UCIWWEIHR FITTING MODEL","text":"This package has a way to sample from a posterior or prior that is defined in the future paper using the uciwweihr_fit.jl and uciwweihr_model.jl.  We can then generate desired quantities and forecast for a given time period with the posterior predictive distribution, using uciwweihr_gq_pp.jl.  We first generate data using the generate_simulation_data_uciwweihr function which is a non-mispecified version of the model.","category":"page"},{"location":"tutorials/uciwweihr_model_fitting/#1.-Data-Generation.","page":"UCIWWEIHR FITTING MODEL","title":"1. Data Generation.","text":"","category":"section"},{"location":"tutorials/uciwweihr_model_fitting/","page":"UCIWWEIHR FITTING MODEL","title":"UCIWWEIHR FITTING MODEL","text":"using UCIWWEIHR\n# Running simulation function with defaults\ndf = generate_simulation_data_uciwweihr()\nfirst(df, 5)","category":"page"},{"location":"tutorials/uciwweihr_model_fitting/#2.-Sampling-from-the-Posterior-Distribution-and-Posterior-Predictive-Distribution.","page":"UCIWWEIHR FITTING MODEL","title":"2. Sampling from the Posterior Distribution and Posterior Predictive Distribution.","text":"","category":"section"},{"location":"tutorials/uciwweihr_model_fitting/","page":"UCIWWEIHR FITTING MODEL","title":"UCIWWEIHR FITTING MODEL","text":"Here we sample from the posterior distribution using the uciwweihr_fit.jl function.  First, we setup some presets, then have an array where index 1 contains the posterior/prior predictive samples, index 2 contains the posterior/prior generated quantities samples, and index 3 contains the original sampled parameters for the model.","category":"page"},{"location":"tutorials/uciwweihr_model_fitting/","page":"UCIWWEIHR FITTING MODEL","title":"UCIWWEIHR FITTING MODEL","text":"data_hosp = df.hosp\ndata_wastewater = df.log_ww_conc\nobstimes = df.obstimes\nparam_change_times = 1:7:length(obstimes) # Change every week\npriors_only = false\nn_samples = 50\n\nsamples = uciwweihr_fit(\n    data_hosp,\n    data_wastewater,\n    obstimes,\n    param_change_times,\n    priors_only,\n    n_samples\n)\nmodel_output = uciwweihr_gq_pp(\n    samples,\n    data_hosp,\n    data_wastewater,\n    obstimes,\n    param_change_times\n)\n\nfirst(model_output[1][:,1:5], 5)","category":"page"},{"location":"tutorials/uciwweihr_model_fitting/","page":"UCIWWEIHR FITTING MODEL","title":"UCIWWEIHR FITTING MODEL","text":"first(model_output[2][:,1:5], 5)","category":"page"},{"location":"tutorials/uciwweihr_model_fitting/","page":"UCIWWEIHR FITTING MODEL","title":"UCIWWEIHR FITTING MODEL","text":"first(model_output[3][:,1:5], 5)","category":"page"},{"location":"tutorials/uciwweihr_model_fitting/#3.-MCMC-Diagnostic-Plots/Results-Along-with-Posterior-Predictive-Distribution.","page":"UCIWWEIHR FITTING MODEL","title":"3. MCMC Diagnostic Plots/Results Along with Posterior Predictive Distribution.","text":"","category":"section"},{"location":"tutorials/uciwweihr_model_fitting/","page":"UCIWWEIHR FITTING MODEL","title":"UCIWWEIHR FITTING MODEL","text":"We also provide a very basic way to visualize some MCMC diagnostics along with effective sample sizes of desired generated quantities(does not include functionality for time-varying quantities).  Along with this, we can also visualize the posterior predictive distribution with actual observed values, which can be used to examine forecasts generated by the model.","category":"page"},{"location":"tutorials/uciwweihr_model_fitting/","page":"UCIWWEIHR FITTING MODEL","title":"UCIWWEIHR FITTING MODEL","text":"uciwweihr_visualizer(gq_samples = model_output[2], \n                        actual_rt_vals = df.rt, \n                        actual_w_t = df.wt, \n                        save_plots = true)","category":"page"},{"location":"tutorials/uciwweihr_model_fitting/#3.1.-MCMC-Diagnostic-Plots.","page":"UCIWWEIHR FITTING MODEL","title":"3.1. MCMC Diagnostic Plots.","text":"","category":"section"},{"location":"tutorials/uciwweihr_model_fitting/","page":"UCIWWEIHR FITTING MODEL","title":"UCIWWEIHR FITTING MODEL","text":"(Image: Plot 1)","category":"page"},{"location":"tutorials/uciwweihr_model_fitting/#3.2.-Time-Varying-Parameter-Results-Plot.","page":"UCIWWEIHR FITTING MODEL","title":"3.2. Time Varying Parameter Results Plot.","text":"","category":"section"},{"location":"tutorials/uciwweihr_model_fitting/","page":"UCIWWEIHR FITTING MODEL","title":"UCIWWEIHR FITTING MODEL","text":"(Image: Plot 2)","category":"page"},{"location":"tutorials/uciwweihr_model_fitting/#[Tutorial-Contents](@ref-tutorial_home)","page":"UCIWWEIHR FITTING MODEL","title":"Tutorial Contents","text":"","category":"section"},{"location":"tutorials/getting_started/#getting_started","page":"GETTING STARTED","title":"Getting Started","text":"","category":"section"},{"location":"tutorials/getting_started/","page":"GETTING STARTED","title":"GETTING STARTED","text":"Welcome to the Tutorials page for the UCIWWEIHR.jl project. This section provides step-by-step guides and examples to help you get started with using the package and understanding its features.","category":"page"},{"location":"tutorials/getting_started/#1.-Installation.","page":"GETTING STARTED","title":"1. Installation.","text":"","category":"section"},{"location":"tutorials/getting_started/","page":"GETTING STARTED","title":"GETTING STARTED","text":"To install the UCIWWEIHR.jl package, open the Julia REPL and run the following command:","category":"page"},{"location":"tutorials/getting_started/","page":"GETTING STARTED","title":"GETTING STARTED","text":"using Pkg\nPkg.add(\"UCIWWEIHR\")","category":"page"},{"location":"tutorials/getting_started/#[Tutorial-Contents](@ref-tutorial_home)","page":"GETTING STARTED","title":"Tutorial Contents","text":"","category":"section"},{"location":"news/#package-development-news","page":"NEWS","title":"Package Development News","text":"","category":"section"},{"location":"news/#Latest-Updates","page":"NEWS","title":"Latest Updates","text":"","category":"section"},{"location":"news/#Upcoming-Features","page":"NEWS","title":"Upcoming Features","text":"","category":"section"},{"location":"news/#Bug-Fixes-in-Progress","page":"NEWS","title":"Bug Fixes in Progress","text":"","category":"section"},{"location":"news/#Release-Notes","page":"NEWS","title":"Release Notes","text":"","category":"section"},{"location":"reference/#reference","page":"REFERENCE","title":"Reference","text":"","category":"section"},{"location":"reference/#Contents","page":"REFERENCE","title":"Contents","text":"","category":"section"},{"location":"reference/","page":"REFERENCE","title":"REFERENCE","text":"Pages = [\"reference.md\"]","category":"page"},{"location":"reference/#Index","page":"REFERENCE","title":"Index","text":"","category":"section"},{"location":"reference/","page":"REFERENCE","title":"REFERENCE","text":"Pages = [\"reference.md\"]","category":"page"},{"location":"reference/","page":"REFERENCE","title":"REFERENCE","text":"Modules = [UCIWWEIHR]","category":"page"},{"location":"reference/#UCIWWEIHR.ChainsCustomIndex-Tuple{MCMCChains.Chains, BitMatrix}","page":"REFERENCE","title":"UCIWWEIHR.ChainsCustomIndex","text":"ChainsCustomIndex(c::Chains, indices_to_keep::BitMatrix)\n\nReduce Chains object to only wanted indices. \n\nFunction created by Damon Bayer. \n\n\n\n\n\n","category":"method"},{"location":"reference/#UCIWWEIHR.NegativeBinomial2-Tuple{Any, Any}","page":"REFERENCE","title":"UCIWWEIHR.NegativeBinomial2","text":"Create a re-parametrized negative binomial distribution in terms of mean and overdispersion.\n\nArguments\n\nμ: Mean of the distribution.\nϕ: Overdispersion parameter.\n\nReturns\n\nA Distributions.NegativeBinomial distribution object.\n\n\n\n\n\n","category":"method"},{"location":"reference/#UCIWWEIHR.calculate_quantiles-NTuple{4, Any}","page":"REFERENCE","title":"UCIWWEIHR.calculate_quantiles","text":"calculate_quantiles(df, chain, var_prefix, quantiles)\n\nCalculate quantiles for a given chain and variable prefix.  Quantiles can be any user desired quantile.\n\nFunction created by Christian Bernal Zelaya. \n\n\n\n\n\n","category":"method"},{"location":"reference/#UCIWWEIHR.eihr_ode!-NTuple{4, Any}","page":"REFERENCE","title":"UCIWWEIHR.eihr_ode!","text":"eihr_ode!(du, u, p, t)\n\nCalculate the ordinary differential equations (ODEs) for the EIHR model.\n\nParameters:\n\ndu: Array{Float64,1} - The derivative of the state variables.\nu: Array{Float64,1} - The current state variables.\np: Tuple{Float64,Float64,Float64,Float64,Float64} - The model parameters (alpha, gamma, nu, w, epsilon).\nt: Float64 - The current time.\n\nReturns:\n\ndu: Array{Float64,1} - The derivative of the state variables.\n\n\n\n\n\n","category":"method"},{"location":"reference/#UCIWWEIHR.generate_colors-Tuple{Any}","page":"REFERENCE","title":"UCIWWEIHR.generate_colors","text":"generate_ribbon_colors(number_of_colors)\n\nGenerates a vector with colors for ribbons in plots.\n\nFunction created by Christian Bernal Zelaya. \n\n\n\n\n\n","category":"method"},{"location":"reference/#UCIWWEIHR.generate_simulation_data_agent","page":"REFERENCE","title":"UCIWWEIHR.generate_simulation_data_agent","text":"Generating Simulation Data for Agent Based Model\n\nTo generate simulation data using the agent based model, you can use the generate_simulation_data_agent function defined in the UCIWWEIHR.jl package. This function allows you to customize various parameters for the simulation. NOT FINISHED, STILL NEEDS WW AND RT\n\nFunction Signature\n\nArguments\n\nseed::Int64: Seed for random number generation. Default value is 1.\npop_size::Int64: Size of the population. Default value is 1000.\nI_init::Int64: Initial number of infected individuals. Default value is 200.\nH_init::Int64: Initial number of hospitalized individuals. Default value is 20.\nbeta::Float64: Transmission rate. Default value is 0.001.\ngamma::Float64: Rate of exposed individuals becoming infectious. Default value is 1/4.\nnu::Float64: Rate of infected individuals recovering or getting hospitalized. Default value is 1/7.\nepsilon::Float64: Rate of hospitalized individuals recovering. Default value is 1/5.\nw_init::Float64: Probability of an infected individual becoming hospitalized. Default value is 0.35.\n\nReturns\n\ndf::DataFrame: A DataFrame containing the simulation data with columns Time, S, E, I, H, and R.\n\n\n\n\n\n","category":"function"},{"location":"reference/#UCIWWEIHR.generate_simulation_data_uciwweihr","page":"REFERENCE","title":"UCIWWEIHR.generate_simulation_data_uciwweihr","text":"Generating Simulation Data for UCIWWEIHR ODE Compartmental Based Model\n\nTo generate simulation data using the UCIWWEIHR ODE compartmental based model, you can use the generate_simulation_data_uciwweihr function defined in the UCIWWEIHR.jl package. This function allows you to customize various parameters for the simulation.\n\nFunction Signature\n\nArguments\n\ntime_points::Int64: Number of time points wanted for simulation. Default value is 150.\nseed::Int64: Seed for random number generation. Default value is 1.\nE_init::Int64: Initial number of exposed individuals. Default value is 200.\nI_init::Int64: Initial number of infected individuals. Default value is 100.\nH_init::Int64: Initial number of hospitalized individuals. Default value is 20.\ngamma::Float64: Rate of incubation. Default value is 1/4.\nnu::Float64: Rate of leaving the infected compartment. Default value is 1/7.\nepsilon::Float64: Rate of hospitalization recovery. Default value is 1/5.\nrho_gene::Float64: Contribution of infected individual's pathogen genome into wastewater. Default value is 0.011.\ntau::Float64: Scale/variation of the log concentration of pathogen genome in wastewater. Default value is 0.1.\ndf::Float64: Degrees of freedom for generalized t distribution for log concentration of pathogen genome in wastewater. Default value is 29.\nsigma_hosp::Float64: Standard deviation for the negative binomial distribution for hospital data. Default value is 800.\nRt_init::Float64: Initial value of the time-varying reproduction number. Default value is 1.\nsigma_Rt::Float64: Standard deviation for random walk of time-varying reproduction number. Default value is sqrt(0.02).\nw_init::Float64: Initial value of the time-varying hospitalization rate. Default value is 0.35.\nsigma_w::Float64: Standard deviation for random walk of time-varying hospitalization rate. Default value is sqrt(0.02).\n\nReturns\n\ndf::DataFrame: A DataFrame containing the simulation data with columns obstimes, log_ww_conc, hosp, and rt.\n\n\n\n\n\n","category":"function"},{"location":"reference/#UCIWWEIHR.power-Tuple{Any, Any}","page":"REFERENCE","title":"UCIWWEIHR.power","text":"power(a,b)\n\nRaise `a` to the `b` power\n\n\n\n\n\n","category":"method"},{"location":"reference/#UCIWWEIHR.save_plots_to_docs-Tuple{Any, Any}","page":"REFERENCE","title":"UCIWWEIHR.save_plots_to_docs","text":"save_plots_to_docs(plot, filename; format = \"png\")\n\nSaves plots to docs/plots directory.\n\nFunction created by Christian Bernal Zelaya. \n\n\n\n\n\n","category":"method"},{"location":"reference/#UCIWWEIHR.startswith_any-Tuple{Any, Any}","page":"REFERENCE","title":"UCIWWEIHR.startswith_any","text":"startswith_any(name, patterns)\n\nChecks if the name of time varying paramter starts with any of the patterns.\n\nFunction created by Christian Bernal Zelaya. \n\n\n\n\n\n","category":"method"},{"location":"reference/#UCIWWEIHR.uciwweihr_fit","page":"REFERENCE","title":"UCIWWEIHR.uciwweihr_fit","text":"uciwweihr_fit(...)\n\nThis is the sampler for the bayesian semi-parametric model for the wastewater EIHR compartmental model.   The defaults for this fuction will follow those of the default simulation in generatesimulationdatawweihr.jl function.\n\nArguments\n\ndata_hosp: An array of hospital data.\ndata_wastewater: An array of pathogen genome concentration in localized wastewater data.\nobstimes: An array of timepoints for observed hosp/wastewater.\npriors_only::Bool=false: A boolean to indicate if only priors are to be sampled.\nn_samples::Int64=500: Number of samples to be drawn.\nn_chains::Int64=1: Number of chains to be run.\nseed::Int64=2024: Seed for the random number generator.\nE_init_sd::Float64=50.0: Standard deviation for the initial number of exposed individuals.\nE_init_mean::Int64=200: Mean for the initial number of exposed individuals.\nI_init_sd::Float64=20.0: Standard deviation for the initial number of infected individuals.\nI_init_mean::Int64=100: Mean for the initial number of infected individuals.\nH_init_sd::Float64=5.0: Standard deviation for the initial number of hospitalized individuals.\nH_init_mean::Int64=20: Mean for the initial number of hospitalized individuals.\ngamma_sd::Float64=0.02: Standard deviation for the rate of incubation.\nlog_gamma_mean::Float64=log(1/4): Mean for the rate of incubation on log scale.\nnu_sd::Float64=0.02: Standard deviation for the rate of leaving the infected compartment.\nlog_nu_mean::Float64=log(1/7): Mean for the rate of leaving the infected compartment on the log scale.\nepsilon_sd::Float64=0.02: Standard deviation for the rate of hospitalization recovery.\nlog_epsilon_mean::Float64=log(1/5): Mean for the rate of hospitalization recovery on the log scale.\nrho_gene_sd::Float64=0.02: Standard deviation for the rho prior.\nlog_rho_gene_mean::Float64=log(0.011): Mean for the row prior on log scale.\ntau_sd::Float64=0.02: Standard deviation for the scale/variation of the log scale data.\nlog_tau_mean::Float64=log(0.1): Mean for the scale/variation of the log scale data on log scale itself.\ndf_shape::Float64=2.0: Shape parameter for the gamma distribution.\ndf_scale::Float64=10.0: Scale parameter for the gamma distribution.\nsigma_hosp_sd::Float64=50.0: Standard deviation for the negative binomial distribution for hospital data.\nsigma_hosp_mean::Float64=500.0: Mean for the negative binomial distribution for hospital data.\nRt_init_sd::Float64=0.3: Standard deviation for the initial value of the time-varying reproduction number.\nRt_init_mean::Float64=0.2: Mean for the initial value of the time-varying reproduction number.\nsigma_Rt_sd::Float64=0.2: Standard deviation for normal prior of log time-varying reproduction number standard deviation.\nsigma_Rt_mean::Float64=-3.0: Mean for normal prior of log time-varying reproduction number standard deviation.\nw_init_sd::Float64=0.1: Standard deviation for the initial value of the time-varying hospitalization rate.\nw_init_mean::Float64=log(0.35): Mean for the initial value of the time-varying hospitalization rate.\nsigma_w_sd::Float64=0.2: Standard deviation for normal prior of log time-varying hospitalization rate standard deviation.\nsigma_w_mean::Float64=-3.5: Mean for normal prior of time-varying hospitalization rate standard deviation.\nparam_change_times::Array{Float64}: An array of timepoints where the parameters change.\n\nReturns\n\nSamples from the posterior or prior distribution.\n\n\n\n\n\n","category":"function"},{"location":"reference/#UCIWWEIHR.uciwweihr_model","page":"REFERENCE","title":"UCIWWEIHR.uciwweihr_model","text":"uciwweihr_model(...)\n\nThis is the bayesian semi-parametric model for the wastewater EIHR compartmental model.   The defaults for this fuction will follow those of the default simulation in generatesimulationdatawweihr.jl function.\n\nArguments\n\ndata_hosp: An array of hospital data.\ndata_wastewater: An array of pathogen genome concentration in localized wastewater data.\nobstimes: An array of timepoints for observed hosp/wastewater.\nparam_change_times: An array of timepoints where the parameters change.\nE_init_sd::Float64=50.0: Standard deviation for the initial number of exposed individuals.\nE_init_mean::Int64=200: Mean for the initial number of exposed individuals.\nI_init_sd::Float64=20.0: Standard deviation for the initial number of infected individuals.\nI_init_mean::Int64=100: Mean for the initial number of infected individuals.\nH_init_sd::Float64=5.0: Standard deviation for the initial number of hospitalized individuals.\nH_init_mean::Int64=20: Mean for the initial number of hospitalized individuals.\ngamma_sd::Float64=0.02: Standard deviation for the rate of incubation.\nlog_gamma_mean::Float64=log(1/4): Mean for the rate of incubation on log scale.\nnu_sd::Float64=0.02: Standard deviation for the rate of leaving the infected compartment.\nlog_nu_mean::Float64=log(1/7): Mean for the rate of leaving the infected compartment on the log scale.\nepsilon_sd::Float64=0.02: Standard deviation for the rate of hospitalization recovery.\nlog_epsilon_mean::Float64=log(1/5): Mean for the rate of hospitalization recovery on the log scale.\nrho_gene_sd::Float64=0.02: Standard deviation for the rho prior.\nlog_rho_gene_mean::Float64=log(0.011): Mean for the row prior on log scale.\ntau_sd::Float64=0.02: Standard deviation for the scale/variation of the log scale data.\nlog_tau_mean::Float64=log(0.1): Mean for the scale/variation of the log scale data on log scale itself.\ndf_shape::Float64=2.0: Shape parameter for the gamma distribution.\ndf_scale::Float64=10.0: Scale parameter for the gamma distribution.\nsigma_hosp_sd::Float64=50.0: Standard deviation for the negative binomial distribution for hospital data.\nsigma_hosp_mean::Float64=500.0: Mean for the negative binomial distribution for hospital data.\nRt_init_sd::Float64=0.3: Standard deviation for the initial value of the time-varying reproduction number.\nRt_init_mean::Float64=0.2: Mean for the initial value of the time-varying reproduction number.\nsigma_Rt_sd::Float64=0.2: Standard deviation for normal prior of log time-varying reproduction number standard deviation.\nsigma_Rt_mean::Float64=-3.0: Mean for normal prior of log time-varying reproduction number standard deviation.\nw_init_sd::Float64=0.1: Standard deviation for the initial value of the time-varying hospitalization rate.\nw_init_mean::Float64=log(0.35): Mean for the initial value of the time-varying hospitalization rate.\nsigma_w_sd::Float64=0.2: Standard deviation for normal prior of log time-varying hospitalization rate standard deviation.\nsigma_w_mean::Float64=-3.5: Mean for normal prior of time-varying hospitalization rate standard deviation.\n\n\n\n\n\n","category":"function"},{"location":"reference/#UCIWWEIHR.uciwweihr_visualizer-Tuple{}","page":"REFERENCE","title":"UCIWWEIHR.uciwweihr_visualizer","text":"uciwweihr_visualizer(...)\n\nDefault visualizer for results of the UCIWWEIHR model, includes posterior/priors of generated quantities and posterior predictive samples for forecasting.  Forecasting plots will have the observed data alongside.\n\nArguments\n\npp_sampeles: Posterior predictive samples from the posterior/prior distribution, index 1 in uciwweihrgqpp output.\ngq_samples: Generated quantities samples from the posterior/prior distribution, index 2 in uciwweihrgqpp output.\ndata_hosp: An array of hospital data.\ndata_wastewater: An array of pathogen genome concentration in localized wastewater data.\nobstimes: An array of timepoints for observed hosp/wastewater.\nparam_change_times: An array of timepoints where the parameters change.\nactual_rt_vals: An array of actual Rt values if user has access to them assumed to be on a daily scale.  This typically will come from some simulation.  Default is nothing.\nactual_w_t: An array of actual w_t values if user has access to them assumed to be on a daily scale.  This typically will come from some simulation.  Default is nothing.\ndesired_params: A list of lists of parameters to visualize. Each list will be visualized in a separate plot. Default is [[\"Einit\", \"Iinit\", \"Hinit\"], [\"gamma\", \"nu\", \"epsilon\"], [\"rhogene\", \"tau\", \"df\"], [\"sigma_hosp\"]].\ntime_varying_params: A list of time varying parameters to visualize. Default is [\"rtvals\", \"wt\"].\nquantiles: A list of quantiles to calculate for ploting uncertainty. Default is [0.5, 0.8, 0.95].\nmcmcdaigs::Bool=true: A boolean to indicate if user wants to visualize mcmc diagnosis plots and Effective Sample Size(ESS).\ntime_varying_plots::Bool=true: A boolean to indicate if user wants to visualize time varying parameters.    \nsave_plots::Bool=false: A boolean to indicate if user wants to save the plots as pngs into a plots folder.\n\n\n\n\n\n","category":"method"},{"location":"#UCIWWEIHR.jl","page":"HOME","title":"UCIWWEIHR.jl","text":"","category":"section"},{"location":"","page":"HOME","title":"HOME","text":"Welcome to the UCIWWEIHR.jl package documentation!  CURRENTLY UNDER DEVELOPMENT!!!","category":"page"},{"location":"","page":"HOME","title":"HOME","text":"By : Christian O. Bernal Zelaya & Volodymyr M. Minin.","category":"page"},{"location":"#Introduction","page":"HOME","title":"Introduction","text":"","category":"section"},{"location":"","page":"HOME","title":"HOME","text":"UCIWWEIHR.jl is a Julia package for forecasting and nowcasting COVID-19 hospitalizations using pathogen genome concentrations. UCIWWEIHR.jl Package License","category":"page"},{"location":"#Features","page":"HOME","title":"Features","text":"","category":"section"},{"location":"","page":"HOME","title":"HOME","text":"Comprehensive library of Bayesian models\nSimulation of infectious disease \nParameter estimation and model calibration\nSensitivity analysis and uncertainty quantification?\nVisualization of results and MCMC evaluations(live?)","category":"page"},{"location":"#Installation","page":"HOME","title":"Installation","text":"","category":"section"},{"location":"","page":"HOME","title":"HOME","text":"To install UCIWWEIHR.jl, you can use the Julia package manager. Open the Julia REPL and run the following command:","category":"page"},{"location":"","page":"HOME","title":"HOME","text":"using Pkg\nPkg.add(\"UCIWWEIHR\")","category":"page"},{"location":"tutorials/uciwweihr_simulation_data/#uciwweihr_simulation_data","page":"UCIWWEIHR SIMULATION DATA","title":"Generating simulated data with UCIWWEIHR ODE compartmental based model.","text":"","category":"section"},{"location":"tutorials/uciwweihr_simulation_data/","page":"UCIWWEIHR SIMULATION DATA","title":"UCIWWEIHR SIMULATION DATA","text":"This package provides a way to also simulate data using the UCIWWEIHR ODE compartmental based model specified in the future paper.  The function called generate_simulation_data_uciwweihr.jl can be used to generate synthetic data for a given number of samples and features.  Here we provide a demonstration using the default settings of generate_simulation_data_uciwweihr.jl :","category":"page"},{"location":"tutorials/uciwweihr_simulation_data/#1.-Functionality.","page":"UCIWWEIHR SIMULATION DATA","title":"1. Functionality.","text":"","category":"section"},{"location":"tutorials/uciwweihr_simulation_data/","page":"UCIWWEIHR SIMULATION DATA","title":"UCIWWEIHR SIMULATION DATA","text":"using UCIWWEIHR\nusing Plots\n# Running simulation function with defaults\ndf = generate_simulation_data_uciwweihr()\nfirst(df, 5)","category":"page"},{"location":"tutorials/uciwweihr_simulation_data/#2.-Visualizing-UCIWWEIHR-model-results.","page":"UCIWWEIHR SIMULATION DATA","title":"2. Visualizing UCIWWEIHR model results.","text":"","category":"section"},{"location":"tutorials/uciwweihr_simulation_data/","page":"UCIWWEIHR SIMULATION DATA","title":"UCIWWEIHR SIMULATION DATA","text":"Here we can make simple plots to visualize the data generated using the Plots  package.","category":"page"},{"location":"tutorials/uciwweihr_simulation_data/#2.1.-Concentration-of-pathogen-genome-in-wastewater(WW).","page":"UCIWWEIHR SIMULATION DATA","title":"2.1. Concentration of pathogen genome in wastewater(WW).","text":"","category":"section"},{"location":"tutorials/uciwweihr_simulation_data/","page":"UCIWWEIHR SIMULATION DATA","title":"UCIWWEIHR SIMULATION DATA","text":"plot(df.obstimes, df.log_ww_conc,\n    label=nothing,\n    xlabel=\"Obstimes\", \n    ylabel=\"Conc. of Pathogen Genome in WW\", \n    title=\"Plot of Conc. of Pathogen Genome in WW Over Time\")","category":"page"},{"location":"tutorials/uciwweihr_simulation_data/#2.2.-Hospitalizations.","page":"UCIWWEIHR SIMULATION DATA","title":"2.2. Hospitalizations.","text":"","category":"section"},{"location":"tutorials/uciwweihr_simulation_data/","page":"UCIWWEIHR SIMULATION DATA","title":"UCIWWEIHR SIMULATION DATA","text":"plot(df.obstimes, df.hosp, \n    label=nothing,\n    xlabel=\"Obstimes\", \n    ylabel=\"Hosp\", \n    title=\"Plot of Hosp Over Time\")","category":"page"},{"location":"tutorials/uciwweihr_simulation_data/#2.3.-Reproductive-number.","page":"UCIWWEIHR SIMULATION DATA","title":"2.3. Reproductive number.","text":"","category":"section"},{"location":"tutorials/uciwweihr_simulation_data/","page":"UCIWWEIHR SIMULATION DATA","title":"UCIWWEIHR SIMULATION DATA","text":"plot(df.obstimes, df.rt, \n    label=nothing,\n    xlabel=\"Obstimes\", \n    ylabel=\"Rt\", \n    title=\"Plot of Rt Over Time\")","category":"page"},{"location":"tutorials/uciwweihr_simulation_data/#2.4.-Hospitalization-rate.","page":"UCIWWEIHR SIMULATION DATA","title":"2.4. Hospitalization rate.","text":"","category":"section"},{"location":"tutorials/uciwweihr_simulation_data/","page":"UCIWWEIHR SIMULATION DATA","title":"UCIWWEIHR SIMULATION DATA","text":"plot(df.obstimes, df.wt, \n    label=nothing,\n    xlabel=\"Obstimes\", \n    ylabel=\"Rt\", \n    title=\"Plot of Hospitalization Rate Over Time\")","category":"page"},{"location":"tutorials/uciwweihr_simulation_data/#[Tutorial-Contents](@ref-tutorial_home)","page":"UCIWWEIHR SIMULATION DATA","title":"Tutorial Contents","text":"","category":"section"},{"location":"package_development/#package-development-updates","page":"PACKAGE DEVELOPMENT","title":"Package Development Updates","text":"","category":"section"},{"location":"package_development/#[Version-X.X.X]-YYYY-MM-DD","page":"PACKAGE DEVELOPMENT","title":"[Version X.X.X] - YYYY-MM-DD","text":"","category":"section"},{"location":"package_development/#Added","page":"PACKAGE DEVELOPMENT","title":"Added","text":"","category":"section"},{"location":"package_development/","page":"PACKAGE DEVELOPMENT","title":"PACKAGE DEVELOPMENT","text":"Feature 1: Description of the feature that was added.\nFeature 2: Description of the feature that was added.","category":"page"},{"location":"package_development/#Changed","page":"PACKAGE DEVELOPMENT","title":"Changed","text":"","category":"section"},{"location":"package_development/","page":"PACKAGE DEVELOPMENT","title":"PACKAGE DEVELOPMENT","text":"Change 1: Description of the change that was made.\nChange 2: Description of the change that was made.","category":"page"},{"location":"package_development/#Deprecated","page":"PACKAGE DEVELOPMENT","title":"Deprecated","text":"","category":"section"},{"location":"package_development/","page":"PACKAGE DEVELOPMENT","title":"PACKAGE DEVELOPMENT","text":"Deprecated feature 1: Description of the feature that was deprecated.\nDeprecated feature 2: Description of the feature that was deprecated.","category":"page"},{"location":"package_development/#Removed","page":"PACKAGE DEVELOPMENT","title":"Removed","text":"","category":"section"},{"location":"package_development/","page":"PACKAGE DEVELOPMENT","title":"PACKAGE DEVELOPMENT","text":"Removed feature 1: Description of the feature that was removed.\nRemoved feature 2: Description of the feature that was removed.","category":"page"},{"location":"package_development/#Fixed","page":"PACKAGE DEVELOPMENT","title":"Fixed","text":"","category":"section"},{"location":"package_development/","page":"PACKAGE DEVELOPMENT","title":"PACKAGE DEVELOPMENT","text":"Bug fix 1: Description of the bug that was fixed.\nBug fix 2: Description of the bug that was fixed.","category":"page"}]
}
