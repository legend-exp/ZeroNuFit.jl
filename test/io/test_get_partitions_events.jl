using Pkg
Pkg.activate(".") # activate the environment
Pkg.instantiate() # instantiate the environment
include("../../src/ZeroNuFit.jl")
using .ZeroNuFit
using TypedTables

@testset "test_get_partitions_events" begin

    @info "Testing function to retrieve input parameters (function 'get_partitions_events' in src/utils.jl)"
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

    part_event_index = nothing
    events = nothing
    partitions = nothing
    fit_ranges = nothing
    try
        part_event_index, events, partitions, fit_ranges =
            ZeroNuFit.Utils.get_partitions_events(config)
    catch e
        @error "Error in get_partitions_events: $e"
        throw(e)
    end

    expected_part_event_index = [1]
    expected_partitions = Table(
        experiment = Array(["GERDA"]),
        fit_group = Array(["all_phase_II"]),
        bkg_name = Array([:B_gerda_all_pII]),
        signal_name = Array([:gaussian]),
        energy_reso_name = Array([:αr_all]),
        energy_bias_name = Array([:αb_all]),
        eff_name = Array([:αe_all]),
        detector = Array(["ANG4"]),
        part_name = Array(["part00"]),
        start_ts = Array([1450622572]),
        end_ts = Array([1469119346]),
        eff_tot = Array([0.476981]),
        eff_tot_sigma = Array([0.0395483]),
        width = Array([1.328561380042463]),
        width_sigma = Array([0.08067940552016985]),
        exposure = Array([0.987039]),
        bias = Array([-0.33885135]),
        bias_sigma = Array([0.07638651]),
        frac = Array([nothing]),
        tau = Array([nothing]),
        sigma = Array([nothing]),
    )
    expected_events = [[1995.2452]]
    expected_fit_ranges = Dict(
        "all_phase_II" =>
            [[1930.0, 2098.511], [2108.511, 2113.513], [2123.513, 2190.0]],
    )


    @testset "Check part_event_index validity" begin
        @test part_event_index == expected_part_event_index
    end
    @testset "Check partitions validity" begin
        @test partitions == expected_partitions
    end
    @testset "Check events validity" begin
        @test events == expected_events
    end
    @testset "Check fit_ranges validity" begin
        @test fit_ranges == expected_fit_ranges
    end

end
