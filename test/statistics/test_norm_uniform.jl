using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("../../src/ZeroNuFit.jl")
using .ZeroNuFit
using Test

@testset "test_norm_uniform" begin

    @info "Testing normalized uniform PDF (function 'norm_uniform' in src/likelihood.jl)"

    # Test case 1: Simple single range [0, 1]
    @testset "Single range [0, 1]" begin
        fit_range = [[0.0, 1.0]]
        p = (B = 1.0,)  # Dummy parameter
        result = ZeroNuFit.Likelihood.norm_uniform(0.5, p, fit_range)
        expected = 1.0  # Uniform over width 1
        @test result ≈ expected
    end

    # Test case 2: Single range with width 10
    @testset "Single range with width 10" begin
        fit_range = [[0.0, 10.0]]
        p = (B = 1.0,)
        result = ZeroNuFit.Likelihood.norm_uniform(5.0, p, fit_range)
        expected = 1.0 / 10.0  # Uniform over width 10
        @test result ≈ expected
    end

    # Test case 3: Multiple disjoint ranges
    @testset "Multiple disjoint ranges" begin
        fit_range = [[0.0, 10.0], [20.0, 30.0]]
        p = (B = 1.0,)
        result = ZeroNuFit.Likelihood.norm_uniform(5.0, p, fit_range)
        expected = 1.0 / 20.0  # Total width is 10 + 10 = 20
        @test result ≈ expected
        
        # Test at a different point
        result2 = ZeroNuFit.Likelihood.norm_uniform(25.0, p, fit_range)
        @test result2 ≈ expected  # Uniform, so same value everywhere
    end

    # Test case 4: Three ranges
    @testset "Three ranges" begin
        fit_range = [[1930.0, 1950.0], [1970.0, 1990.0], [2000.0, 2050.0]]
        p = (B = 1.0,)
        result = ZeroNuFit.Likelihood.norm_uniform(1940.0, p, fit_range)
        total_width = 20.0 + 20.0 + 50.0  # = 90.0
        expected = 1.0 / total_width
        @test result ≈ expected atol=1e-10
    end

    # Test case 5: Check that result is independent of x (uniformity)
    @testset "Uniformity check" begin
        fit_range = [[0.0, 100.0]]
        p = (B = 1.0,)
        result1 = ZeroNuFit.Likelihood.norm_uniform(10.0, p, fit_range)
        result2 = ZeroNuFit.Likelihood.norm_uniform(50.0, p, fit_range)
        result3 = ZeroNuFit.Likelihood.norm_uniform(90.0, p, fit_range)
        @test result1 ≈ result2 ≈ result3
    end

end
