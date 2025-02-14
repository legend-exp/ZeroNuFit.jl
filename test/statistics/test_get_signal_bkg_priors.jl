using Pkg
Pkg.activate(".")
Pkg.instantiate()
using Random
include("../../src/ZeroNuFit.jl")
using .ZeroNuFit
using IntervalSets
using Distributions

Base.exit(code::Int) = throw(ArgumentError("exit code $code"))

@testset "test_get_signal_bkg_priors" begin

    @info "Testing retrieval of signal & background pdfs (function 'get_signal_bkg_priors' in src/likelihood.jl)"
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

    distrS = nothing
    distrB = nothing
    try
        distrS, distrB = ZeroNuFit.Likelihood.get_signal_bkg_priors(config)
    catch e
        @error "Error in 'get_signal_bkg_priors' evaluation: $e"
        throw(e)
    end

    @testset "Check distrS/distrB is valid" begin
        @test !isnothing(distrS)
        @test !isnothing(distrB)
    end

    expected_distrS = IntervalSets.ClosedInterval(0, 1000)
    expected_distrB = IntervalSets.ClosedInterval(0, 0.1)
    @testset "Check distrS/distrB accuracy" begin
        @test distrS == expected_distrS
        @test distrB == expected_distrB
    end

    # not existing signal prior
    config["signal"]["prior"] = "fancy_prior"
    distrS = nothing
    distrB = nothing
    @test_throws ArgumentError ZeroNuFit.Likelihood.get_signal_bkg_priors(config)

    # not existing background prior
    config["signal"]["prior"] = "uniform"
    config["bkg"]["prior"] = "fancy_prior"
    distrS = nothing
    distrB = nothing
    @test_throws ArgumentError ZeroNuFit.Likelihood.get_signal_bkg_priors(config)

    # sqrt signal prior
    config["signal"]["prior"] = "sqrt"
    config["bkg"]["prior"] = "uniform"
    distrS = nothing
    distrB = nothing
    distrS, distrB = ZeroNuFit.Likelihood.get_signal_bkg_priors(config)
    expected_distrS = IntervalSets.ClosedInterval(0, 1000)
    @test distrS == expected_distrS

    # loguniform signal prior
    config["signal"]["prior"] = "loguniform"
    distrS = nothing
    distrB = nothing
    distrS, distrB = ZeroNuFit.Likelihood.get_signal_bkg_priors(config)
    expected_distrS = LogUniform(0.01, 1000)
    @test distrS == expected_distrS
end
