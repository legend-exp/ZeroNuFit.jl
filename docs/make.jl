using Pkg
Pkg.activate(".")
Pkg.instantiate()

using Documenter
push!(LOAD_PATH, "../src")
using ZeroNuFit

makedocs(
    modules = [ZeroNuFit],
    sitename = "ZeroNuFit.jl",
    authors = "S. Calgaro, T. Dixon",
    format = Documenter.HTML(
        #mathengine = MathJax(),
        # Set size thresholds for HTML generation
        size_threshold = 400 * 1024,       # Hard limit (e.g., 400 KiB)
        size_threshold_warn = 300 * 1024  # Warning limit (e.g., 300 KiB)
    ),
    pages = [
        "Home" => "index.md",
        "First steps" => "installation.md",
        "Configuration file" => "config.md",
        "Partitions and events" => "inputs.md",
        "Generating toys" => "toys.md",
        "Tutorial" => "tutorial.md",
    ]
)
