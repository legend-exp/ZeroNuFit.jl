using Pkg
Pkg.activate(".") # activate the environment
Pkg.instantiate() # instantiate the environment
include("../../src/ZeroNuFit.jl")
using .ZeroNuFit

@testset "test_get_settings" begin

    @info "Testing function to retrieve settings information starting from config (function 'get_settings' in src/utils.jl)"
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

    settings = nothing
    try
        settings = ZeroNuFit.Utils.get_settings(config)
    catch e
        @error "Error in get_settings: $e"
        throw(e)
    end

    @testset "Check settings is valid" begin
        @test !isnothing(settings)
    end

    expected_settings = Dict(
        :eff_correlated => true,
        :eff_fixed => true,
        :energy_bias_fixed => true,
        :energy_bias_correlated => false,
        :energy_res_fixed => true,
        :energy_res_correlated => false,
        :bkg_only => false,
    )
    @testset "Check settings accuracy" begin
        @test settings == expected_settings
    end
end
