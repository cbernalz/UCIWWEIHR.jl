using Documenter, ComModeForecasting

DocMeta.setdocmeta!(ComModeForecasting, :DocTestSetup, :(using ComModeForecasting); recursive = true)

makedocs(;
  modules = [ComModeForecasting],
  doctest = true,
  linkcheck = true, # Rely on Lint.yml/lychee for the links
  authors = "Chrisitan O. Bernal Zelaya <cbernalz@uci.edu>",
  repo = "https://github.com/cbernalz/ComModeForecasting.jl/blob/{commit}{path}#{line}",
  sitename = "ComModeForecasting.jl",
  format = Documenter.HTML(;
    prettyurls = get(ENV, "CI", "false") == "true",
    repolink = "https://cbernalz.github.io/ComModeForecasting.jl",
    size_threshold_warn = nothing,
    size_threshold = nothing # Increase this value
  ),
  pages = [
    "HOME" => "index.md", 
    "TUTORIALS" => [
      "TUTORIAL CONTENTS" => "tutorial_index.md",
      "GETTING STARTED" => "tutorials/getting_started.md",
      "ComModeForecasting SIMULATION DATA" => "tutorials/ComModeForecasting_simulation_data.md",
      "ComModeForecasting FITTING MODEL W/OUT FORECASTING" => "tutorials/ComModeForecasting_model_fitting_no_forecast.md",
      "ComModeForecasting FITTING MODEL W/ FORECASTING" => "tutorials/ComModeForecasting_model_fitting_forecast.md",
    ]
    ,
    "REFERENCE" => "reference.md", 
    "LICENSE" => "license.md",
    ],
)

deploydocs(; repo = "github.com/cbernalz/ComModeForecasting.jl", push_preview = false)