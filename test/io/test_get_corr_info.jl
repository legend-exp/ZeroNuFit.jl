using Pkg
Pkg.activate(".")
Pkg.instantiate()
using Random
include("../../src/ZeroNuFit.jl")
using .ZeroNuFit

@testset "test_get_corr_info" begin

    @info "Testing function to retrieve bkg correlation info (function 'get_corr_info' in src/utils.jl)"

    # no entry for correlated bkg
    config = Dict("bkg" => Dict())

    corr = false
    hier_mode = nothing
    hier_range = nothing
    try
        corr, hier_mode, hier_range = ZeroNuFit.Utils.get_corr_info(config)
    catch e
        @error "Error in 'get_corr_info' evaluation: $e"
        throw(e)
    end

    @test corr == false
    @test hier_mode == nothing
    @test hier_range == nothing

    # not-correlated bkg
    config = Dict("bkg" => Dict("correlated" => Dict("mode" => "none", "range" => "none")))

    corr = false
    hier_mode = nothing
    hier_range = nothing
    try
        corr, hier_mode, hier_range = ZeroNuFit.Utils.get_corr_info(config)
    catch e
        @error "Error in 'get_corr_info' evaluation: $e"
        throw(e)
    end

    @test corr == false
    @test hier_mode == nothing
    @test hier_range == nothing

    # correlated bkg
    config = Dict(
        "bkg" => Dict("correlated" => Dict("mode" => "lognormal", "range" => [0, 0.1])),
    )

    corr = false
    hier_mode = nothing
    hier_range = nothing
    try
        corr, hier_mode, hier_range = ZeroNuFit.Utils.get_corr_info(config)
    catch e
        @error "Error in 'get_corr_info' evaluation: $e"
        throw(e)
    end

    expected_value = true
    @testset "Check corr accuracy" begin
        @test corr == expected_value
    end
    expected_value = config["bkg"]["correlated"]["mode"]
    @testset "Check hier_mode accuracy" begin
        @test hier_mode == expected_value
    end
    expected_value = config["bkg"]["correlated"]["range"]
    @testset "Check hier_range accuracy" begin
        @test hier_range == expected_value
    end

end
