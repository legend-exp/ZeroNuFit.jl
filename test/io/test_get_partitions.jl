using Pkg
Pkg.activate(".") # activate the environment
Pkg.instantiate() # instantiate the environment
include("../../src/ZeroNuFit.jl")
using .ZeroNuFit
using TypedTables

@testset "test_get_partitions" begin

    @info "Testing function to retrieve partitions (function 'get_partitions' in src/utils.jl)"
    present_dir = @__DIR__

    # one partition
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

    partitions = nothing
    fit_ranges = nothing
    try
        partitions, fit_ranges = ZeroNuFit.Utils.get_partitions(config)
    catch e
        @error "Error in get_partitions: $e"
        throw(e)
    end

    @testset "Check partitions is valid" begin
        @test !isnothing(partitions)
    end
    @testset "Check fit_ranges is valid" begin
        @test !isnothing(fit_ranges)
    end

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

    @testset "Check partitions accuracy" begin
        @test partitions == expected_partitions
    end

    expected_fit_ranges = Dict(
        "all_phase_II" =>
            [[1930.0, 2098.511], [2108.511, 2113.513], [2123.513, 2190.0]],
    )
    @testset "Check fit_ranges accuracy" begin
        @test fit_ranges == expected_fit_ranges
    end

    # more partitions
    config["partitions"] = [
        joinpath(present_dir, "../inputs/partitions_test.json"),
        joinpath(present_dir, "../inputs/partitions_test_2.json"),
    ]
    partitions = nothing
    fit_ranges = nothing
    partitions, fit_ranges = ZeroNuFit.Utils.get_partitions(config)
    expected_partitions = Table(
        experiment = Array(["GERDA", "MJD"]),
        fit_group = Array(Any["all_phase_II", "mjd-DS0"]),
        bkg_name = Array(Any[:B_gerda_all_pII, Symbol("mjd-DS0")]),
        signal_name = Array([:gaussian, :gaussian_plus_lowEtail]),
        energy_reso_name = Array([:αr_all, :αr_all]),
        energy_bias_name = Array([:αb_all, :αb_all]),
        eff_name = Array([:αe_all, :αe_all]),
        detector = Array(["ANG4", "DS0"]),
        part_name = Array(["part00", "ds0"]),
        start_ts = Array([1450622572, 0]),
        end_ts = Array([1469119346, 5]),
        eff_tot = Array([0.476981, 0.2]),
        eff_tot_sigma = Array([0.0395483, 0.01]),
        width = Array([1.328561380042463, 1.05]),
        width_sigma = Array([0.08067940552016985, 0.04]),
        exposure = Array([0.987039, 1.08]),
        bias = Array([-0.33885135, 0.02]),
        bias_sigma = Array([0.07638651, 0.08]),
        frac = Array([nothing, 0.1]),
        tau = Array([nothing, 1.05]),
        sigma = Array([nothing, 1.1]),
    )
    expected_fit_ranges = Dict(
        "all_phase_II" =>
            [[1930.0, 2098.511], [2108.511, 2113.513], [2123.513, 2190.0]],
        "mjd-DS0" => [
            [1950.0, 2098.511],
            [2108.511, 2113.513],
            [2123.513, 2199.1],
            [2209.1, 2350.0],
        ],
    )
    #@test partitions == expected_partitions
    @test partitions == expected_partitions
    @test fit_ranges == expected_fit_ranges

end
