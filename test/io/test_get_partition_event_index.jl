using Pkg
Pkg.activate(".") # activate the environment
Pkg.instantiate() # instantiate the environment

include("../../src/ZeroNuFit.jl")
using .ZeroNuFit

@testset "test_get_partition_event_index" begin

    @info "Testing function to retrieve index of events given partitions (function 'get_partition_event_index' in src/utils.jl)"
    present_dir = @__DIR__

    config = Dict(
        "bkg_only" => false,
        "bkg" => Dict(
            "correlated" => Dict("range" => "none", "mode" => "none"),
            "prior" => "uniform",
            "upper_bound" => 0.1,
        ),
        "signal" => Dict("prior" => "uniform", "upper_bound" => 1000),
        "events" => [joinpath(present_dir, "../inputs/events_test.json")],
        "partitions" => [joinpath(present_dir, "../inputs/partitions_test.json")],
        "output_path" => "tests",
        "bat_fit" => Dict("nsteps" => 10000.0, "nchains" => 4),
        "nuisance" => Dict(
            "efficiency" => Dict("fixed" => true, "correlated" => true),
            "energy_bias" => Dict("fixed" => true, "correlated" => false),
            "energy_res" => Dict("fixed" => true, "correlated" => false),
        ),
        "plot" => Dict(
            "bandfit_and_data" => false,
            "alpha" => 0.3,
            "fit_and_data" => false,
            "scheme" => "blue",
        ),
    )

    # 1 event case
    partitions = nothing
    part_event_index = nothing
    fit_ranges = nothing
    partitions, fit_ranges = ZeroNuFit.Utils.get_partitions(config)
    events = ZeroNuFit.Utils.get_events(config["events"][1], partitions)
    try
        part_event_index = ZeroNuFit.Utils.get_partition_event_index(events, partitions)
    catch e
        @error "Error in get_partition_event_index: $e"
        throw(e)
    end

    @testset "Check part_event_index is valid" begin
        @test !isnothing(part_event_index)
    end

    expected_part_event_index = [1]
    @testset "Check part_event_index accuracy (1 event)" begin
        @test part_event_index == expected_part_event_index
    end

    # no events case
    events = Array{Vector{Float64}}(undef, length(partitions))
    for (idx, part) in enumerate(partitions)
        events[idx] = Vector{Float64}[]
    end
    part_event_index = ZeroNuFit.Utils.get_partition_event_index(events, partitions)
    @testset "Check part_event_index accuracy (0 events)" begin
        @test part_event_index == [0]
    end

end
