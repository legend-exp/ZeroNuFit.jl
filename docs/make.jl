using Pkg
Pkg.activate(".")
Pkg.instantiate()

using Documenter
push!(LOAD_PATH, "../src")
using ZeroNuFit

makedocs(
    debug = true,
    modules = [ZeroNuFit],
    sitename = "ZeroNuFit.jl",
    authors = "Sofia Calgaro, Toby Dixon, Luigi Pertoldi",
    format = Documenter.HTML(
        # size thresholds for HTML generation
        size_threshold = 400 * 1024,       # hard limit
        size_threshold_warn = 300 * 1024,  # warning limit 
    ),
    pages = [
        "Home" => "index.md",
        "Likelihood implementation" => "likelihood.md",
        "First steps" => "installation.md",
        "Configuration file" => "config.md",
        "Partitions and events" => "inputs.md",
        "Generating toys" => "toys.md",
        "Tutorial" => "tutorial.md",
        "API" => "api.md",
    ],
)
