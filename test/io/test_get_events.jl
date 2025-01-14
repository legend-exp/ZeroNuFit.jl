using Pkg
Pkg.activate(".") # activate the environment
Pkg.instantiate() # instantiate the environment
include("../../src/ZeroNuFit.jl")
using .ZeroNuFit
include("../../main.jl")
include("../../src/utils.jl")

Base.exit(code::Int) = throw(ArgumentError("exit code $code"))

@testset "test_get_events" begin

    @info "Testing function to retrieve events given partitions (function 'get_events' in src/utils.jl)"
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
            "energy_scale" => Dict("fixed" => true, "correlated" => false),
        ),
        "plot" => Dict(
            "bandfit_and_data" => false,
            "alpha" => 0.3,
            "fit_and_data" => false,
            "scheme" => "blue",
        ),
    )

    partitions = nothing
    events = nothing

    partitions, fit_ranges = ZeroNuFit.get_partitions(config)
    try
        events = ZeroNuFit.get_events(config["events"][1], partitions)
    catch e
        @error "Error in get_events: $e"
        throw(e)
    end

    @testset "Check events is valid" begin
        @test !isnothing(events)
    end

    expected_events = [[1995.2452]]
    @testset "Check events accuracy" begin
        @test events == expected_events
    end

    # event with no partition
    @test_throws ArgumentError ZeroNuFit.get_events(
        joinpath(present_dir, "../inputs/events_fake.json"),
        partitions,
    )

    # not-existing file
    @test_throws ArgumentError ZeroNuFit.get_events("not_existing_file.json", partitions)
end
