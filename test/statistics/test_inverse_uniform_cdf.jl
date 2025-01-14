using Pkg
Pkg.activate(".")
Pkg.instantiate()
using Random
include("../../src/ZeroNuFit.jl")
using .ZeroNuFit
include("../../main.jl")
include("../../src/utils.jl")

@testset "test_inverse_uniform_cdf" begin

    @info "Testing function for getting the inverse CDF (function 'inverse_uniform_cdf' in src/utils.jl)"

    Random.seed!(123)
    n = rand() # 0.521213795535383
    fit_range = [[1, 10], [12, 30]]

    res = nothing
    try
        res = ZeroNuFit.inverse_uniform_cdf(n, fit_range)
    catch e
        @error "Error in 'inverse_uniform_cdf' evaluation: $e"
        throw(e)
    end

    @testset "Check res is valid" begin
        @test !isnothing(res)
    end

    expected_value = 17.072772479455345
    tolerance = 1e-5
    @testset "Check res accuracy" begin
        diff = abs(res - expected_value)
        @test diff <= tolerance
    end

    # lower edge
    res = ZeroNuFit.inverse_uniform_cdf(0, fit_range)
    expected_value = 1
    tolerance = 1e-5
    @testset "Check res accuracy at lower edge" begin
        diff = abs(res - expected_value)
        @test diff <= tolerance
    end
    # upper edge
    res = ZeroNuFit.inverse_uniform_cdf(1, fit_range)
    expected_value = 30
    tolerance = 1e-5
    @testset "Check res accuracy at upper edge" begin
        diff = abs(res - expected_value)
        @test diff <= tolerance
    end

    # simple case
    res = ZeroNuFit.inverse_uniform_cdf(0.5, [[0, 1]])
    expected_value = 0.5
    tolerance = 1e-5
    @testset "Check res accuracy in [0,1] at p=0.5" begin
        diff = abs(res - expected_value)
        @test diff <= tolerance
    end
end
