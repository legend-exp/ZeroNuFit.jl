using Pkg
Pkg.activate(".") # activate the environment
Pkg.instantiate() # instantiate the environment
include("../../src/ZeroNuFit.jl")
using .ZeroNuFit
include("../../main.jl")
include("../../src/plotting.jl")

@testset "test_fit_model" begin

    @info "Testing function to retrieve fit function added to the final plot evaluated at an energy of 2039.05 keV (function 'fit_model' in src/plotting.jl)"

    fit_ranges =
        Dict("l200a" => [[1930.0, 2098.511], [2108.511, 2113.513], [2123.513, 2190.0]])
    part_event_index = [1]
    partitions = Table(
        experiment = Array(["L200"]),
        fit_group = Array(["l200a"]),
        bkg_name = Array([:B_l200a_all]),
        signal_name = Array([:gaussian]),
        energy_reso_name = Array([:αr_all]),
        energy_bias_name = Array([:αb_all]),
        eff_name = Array([:αe_all]),
        detector = Array(["V09374A"]),
        part_name = Array(["part0002"]),
        start_ts = Array([1686522477]),
        end_ts = Array([1690190843]),
        eff_tot = Array([0.5]),
        eff_tot_sigma = Array([0.01]),
        width = Array([1.3]),
        width_sigma = Array([0.01]),
        exposure = Array([1]),
        bias = Array([0.1]),
        bias_sigma = Array([0.01]),
        frac = Array([nothing]),
        tau = Array([nothing]),
        sigma = Array([nothing]),
    )

    settings = Dict()
    settings[:energy_scale_fixed] = true
    settings[:energy_scale_correlated] = false
    settings[:eff_fixed] = true
    settings[:eff_correlated] = false
    settings[:bkg_only] = true

    sqrt_prior = false
    s_max = nothing
    bkg_shape = :uniform
    p = (S = 100, αe_all = 0.1, ω = [1.1], 𝛥 = [0.1], B_l200a_all = 2E-4)
    x = 2039.05

    total = nothing
    try
        total = ZeroNuFit.fit_model(
            part_event_index[1],
            partitions[1],
            p,
            settings,
            bkg_shape,
            fit_ranges[partitions[1].fit_group],
            x,
        )
    catch e
        @error "Error in fit_model: $e"
        throw(e)
    end

    @testset "Check total is valid" begin
        @test !isnothing(total)
    end

    expected_total = 0.0002
    tolerance = 1e-3
    @testset "Check events accuracy" begin
        diff = abs(total - expected_total)
        @test diff <= tolerance
    end

    settings[:bkg_only] = false
    total = nothing
    total = ZeroNuFit.fit_model(
        part_event_index[1],
        partitions[1],
        p,
        settings,
        bkg_shape,
        fit_ranges[partitions[1].fit_group],
        x,
    )
    expected_total = 0.0845283663997091
    tolerance = 1e-3
    @testset "Check events accuracy" begin
        diff = abs(total - expected_total)
        @test diff <= tolerance
    end
end
