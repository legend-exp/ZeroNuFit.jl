
using ZeroNuFit

# List of submodules in ZeroNuFit
submodules = [
    ZeroNuFit.Analysis,
    ZeroNuFit.Constants,
    ZeroNuFit.Fitting,
    ZeroNuFit.Likelihood,
    ZeroNuFit.Plotting,
    ZeroNuFit.Utils,
]

# Function to list functions in a submodule
function get_functions(mod)
    return filter(s -> isa(getfield(mod, s), Function), names(mod))
end

# Collect the functions in each submodule
functions_in_submodules = Dict(mod => get_functions(mod) for mod in submodules)

# Print the functions in each submodule
for (mod, functions) in functions_in_submodules
    println("Functions in $(mod):")
    println(functions)
    println()
end
