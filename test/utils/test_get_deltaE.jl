using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("../../src/ZeroNuFit.jl")
using .ZeroNuFit
using Test

@testset "test_get_deltaE" begin

    @info "Testing function to get net width of fit range (function 'get_deltaE' in src/utils.jl)"

    # Test case 1: Single range
    @testset "Single range" begin
        fit_range = [[1930, 1950]]
        result = ZeroNuFit.Utils.get_deltaE(fit_range)
        expected = 20.0
        @test result ≈ expected
    end

    # Test case 2: Multiple ranges
    @testset "Multiple ranges" begin
        fit_range = [[1930, 1950], [1970, 1990], [2000, 2050]]
        result = ZeroNuFit.Utils.get_deltaE(fit_range)
        expected = 20.0 + 20.0 + 50.0  # = 90.0
        @test result ≈ expected
    end

    # Test case 3: Disjoint ranges with gaps
    @testset "Disjoint ranges" begin
        fit_range = [[0, 10], [20, 30], [50, 60]]
        result = ZeroNuFit.Utils.get_deltaE(fit_range)
        expected = 10.0 + 10.0 + 10.0  # = 30.0
        @test result ≈ expected
    end

    # Test case 4: Single point range (zero width)
    @testset "Zero width range" begin
        fit_range = [[1000, 1000]]
        result = ZeroNuFit.Utils.get_deltaE(fit_range)
        expected = 0.0
        @test result ≈ expected
    end

    # Test case 5: Floating point ranges
    @testset "Floating point ranges" begin
        fit_range = [[1.5, 2.7], [3.1, 4.9]]
        result = ZeroNuFit.Utils.get_deltaE(fit_range)
        expected = 1.2 + 1.8  # = 3.0
        @test result ≈ expected atol = 1e-10
    end

end
