using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("../../src/ZeroNuFit.jl")
using .ZeroNuFit
using Test

@testset "test_get_corr_info" begin

    @info "Testing function to retrieve correlated background info (function 'get_corr_info' in src/utils.jl)"

    # Test case 1: No correlated key in config
    @testset "No correlated key" begin
        config = Dict("bkg" => Dict("upper_bound" => 100))
        corr, hier_mode, hier_range = ZeroNuFit.Utils.get_corr_info(config)
        @test corr == false
        @test hier_mode === nothing
        @test hier_range === nothing
    end

    # Test case 2: Correlated mode is "none"
    @testset "Correlated mode is none" begin
        config =
            Dict("bkg" => Dict("correlated" => Dict("mode" => "none", "range" => [0, 1])))
        corr, hier_mode, hier_range = ZeroNuFit.Utils.get_corr_info(config)
        @test corr == false
        @test hier_mode === nothing
        @test hier_range === nothing
    end

    # Test case 3: Correlated mode is active
    @testset "Correlated mode active" begin
        config = Dict(
            "bkg" => Dict(
                "correlated" => Dict("mode" => "hierarchical", "range" => [0.5, 1.5]),
            ),
        )
        corr, hier_mode, hier_range = ZeroNuFit.Utils.get_corr_info(config)
        @test corr == true
        @test hier_mode == "hierarchical"
        @test hier_range == [0.5, 1.5]
    end

    # Test case 4: Different mode name
    @testset "Different mode name" begin
        config = Dict(
            "bkg" => Dict("correlated" => Dict("mode" => "uniform", "range" => [1, 2])),
        )
        corr, hier_mode, hier_range = ZeroNuFit.Utils.get_corr_info(config)
        @test corr == true
        @test hier_mode == "uniform"
        @test hier_range == [1, 2]
    end

end
