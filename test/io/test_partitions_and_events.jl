using Pkg
Pkg.activate(".") # activate the environment
Pkg.instantiate() # instantiate the environment
using ArgParse
using Logging, LoggingExtras
using JSON
using FilePathsBase
using DensityInterface

include("../../src/ZeroNuFit.jl")
using .ZeroNuFit
include("../../main.jl")
include("../../src/utils.jl")

function main()
    
    config = Dict(
        "bkg_only" => false,
        "bkg" => Dict("correlated" => Dict("range" => "none", "mode" => "none"), "prior" => "uniform", "upper_bound" => 0.1), 
        "signal" => Dict("prior" => "uniform", "upper_bound" => 1000), 
        "events" => ["legend-0vbb-config/tests/events_test.json"], 
        "partitions" => ["legend-0vbb-config/tests/partitions_test.json"], 
        "output_path" => "tests", 
        "bat_fit" => Dict("nsteps" => 10000.0, "nchains" => 4), 
        "nuisance" => Dict("efficiency" => Dict("fixed" => true, "correlated" => true), "energy_scale" => Dict("fixed" => true, "correlated" => false)), 
        "plot" => Dict("bandfit_and_data" => false, "alpha" => 0.3, "fit_and_data" => false, "scheme" => "blue")
    )
    
    part_event_index = nothing
    events = nothing
    partitions = nothing
    fit_ranges = nothing
    try
        part_event_index,events,partitions,fit_ranges= ZeroNuFit.get_partitions_events(config)
        @info "Partitions: $partitions "
        @info "Events: $events"
        @info "Fit ranges: $fit_ranges"
    catch e
        @error "Error in get_partitions_events: $e"
        throw(e)
    end
end
    


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
