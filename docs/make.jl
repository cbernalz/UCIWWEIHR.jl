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
    canonical = "https://cbernalz.github.io/UCIWWEIHR.jl",
  ),
  pages = [
    "HOME" => "index.md", 
    "TUTORIAL" => "tutorial.md",
    "NEWS" => "news.md",
    "PACKAGE DEVELOPMENT" => "package_development.md",
    "REFERENCE" => "reference.md", 
    "LICENSE" => "license.md",
    ],
)

deploydocs(; repo = "github.com/cbernalz/UCIWWEIHR.jl", push_preview = false)