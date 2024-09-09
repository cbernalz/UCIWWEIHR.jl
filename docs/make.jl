using Documenter, UCIWWEIHR

DocMeta.setdocmeta!(UCIWWEIHR, :DocTestSetup, :(using UCIWWEIHR); recursive = true)

makedocs(;
  modules = [UCIWWEIHR],
  doctest = true,
  linkcheck = true, # Rely on Lint.yml/lychee for the links
  authors = "Chrisitan O. Bernal Zelaya <cbernalz@uci.edu>",
  repo = "https://github.com/cbernalz/UCIWWEIHR.jl/blob/{commit}{path}#{line}",
  sitename = "UCIWWEIHR.jl",
  format = Documenter.HTML(;
    prettyurls = get(ENV, "CI", "false") == "true",
    repolink = "https://cbernalz.github.io/UCIWWEIHR.jl",
    size_threshold_warn = nothing,
    size_threshold = nothing # Increase this value
  ),
  pages = [
    "HOME" => "index.md", 
    "TUTORIALS" => [
      "TUTORIAL CONTENTS" => "tutorial_index.md",
      "GETTING STARTED" => "tutorials/getting_started.md",
      "UCIWWEIHR SIMULATION DATA" => "tutorials/uciwweihr_simulation_data.md",
      "AGENT-BASED SIMULATION DATA" => "tutorials/agent_based_simulation_data.md",
      "UCIWWEIHR FITTING MODEL W/OUT FORECASTING" => "tutorials/uciwwiehr_model_fitting_no_forecast.md",
      "UCIWWEIHR FITTING MODEL W/ FORECASTING" => "tutorials/uciwwiehr_model_fitting_forecast.md",
      "UCIWWEIHR REPEATED FORECASTING" => "tutorials/uciwwiehr_model_repeated_forecasts.md",
    ]
    ,
    "REFERENCE" => "reference.md", 
    "LICENSE" => "license.md",
    ],
)

deploydocs(; repo = "github.com/cbernalz/UCIWWEIHR.jl", push_preview = false)