using Pkg
Pkg.activate(".") # activate the environment
Pkg.instantiate() # instantiate the environment

include("../../src/ZeroNuFit.jl")
using .ZeroNuFit
include("../../main.jl")
include("../../src/utils.jl")

@testset "test_get_partition_event_index" begin
    
    @info "Testing function to retrieve index of events given partitions (function 'get_partition_event_index' in src/utils.jl)"
    present_dir = @__DIR__
    
    config = Dict(
        "bkg_only" => false,
        "bkg" => Dict("correlated" => Dict("range" => "none", "mode" => "none"), "prior" => "uniform", "upper_bound" => 0.1), 
        "signal" => Dict("prior" => "uniform", "upper_bound" => 1000), 
        "events" => [joinpath(present_dir, "../inputs/events_test.json")], 
        "partitions" => [joinpath(present_dir, "../inputs/partitions_test.json")], 
        "output_path" => "tests", 
        "bat_fit" => Dict("nsteps" => 10000.0, "nchains" => 4), 
        "nuisance" => Dict("efficiency" => Dict("fixed" => true, "correlated" => true), "energy_scale" => Dict("fixed" => true, "correlated" => false)), 
        "plot" => Dict("bandfit_and_data" => false, "alpha" => 0.3, "fit_and_data" => false, "scheme" => "blue")
    )
    
    partitions = nothing
    part_event_index = nothing
    fit_ranges = nothing
    partitions,fit_ranges = ZeroNuFit.get_partitions(config)
    events = ZeroNuFit.get_events(config["events"][1],partitions)
    try
        part_event_index = ZeroNuFit.get_partition_event_index(events,partitions)
    catch e
        @error "Error in get_partition_event_index: $e"
        throw(e)
    end
    
    @testset "Check part_event_index is valid" begin
        @test !isnothing(part_event_index)
    end
    
    expected_part_event_index=[1]
    @testset "Check part_event_index accuracy" begin
        @test part_event_index == expected_part_event_index
    end
    
end
    
