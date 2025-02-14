using Pkg
Pkg.activate(".")
Pkg.instantiate()
using Random
include("../../src/ZeroNuFit.jl")
using .ZeroNuFit

@testset "test_norm_exponential" begin

    @info "Testing normalised exponential function (function 'norm_exponential' in src/likelihood.jl)"

    # high Rt value
    x = 1965.0
    p = (
        S = 100,
        αe_all = 0.1,
        ω = [1.1],
        𝛥 = [0.1],
        B_l200a_all = 2E-4,
        B_l200a_all_slope = 10,
    )
    b_name = :B_l200a_all
    fit_range = [[1920.0, 1930.0], [1960.0, 1970.0]]
    expo_func = nothing
    try
        expo_func = ZeroNuFit.Likelihood.norm_exponential(x, p, b_name, fit_range)
    catch e
        @error "Error in 'norm_exponential' evaluation: $e"
        throw(e)
    end

    @testset "Check expo_func is valid (high Rt value)" begin
        @test !isnothing(expo_func)
    end

    expected_value = 0.08506327727340002
    tolerance = 1e-3
    @testset "Check expo_func accuracy (high Rt value)" begin
        diff = abs(expo_func - expected_value)
        @test diff <= tolerance
    end

    # low Rt value
    p = (
        S = 100,
        αe_all = 0.1,
        ω = [1.1],
        𝛥 = [0.1],
        B_l200a_all = 2E-4,
        B_l200a_all_slope = 0.00001,
    )
    expo_func = nothing
    try
        expo_func = ZeroNuFit.Likelihood.norm_exponential(x, p, b_name, fit_range)
    catch e
        @error "Error in 'norm_exponential' evaluation: $e"
        throw(e)
    end

    @testset "Check expo_func is valid (low Rt value)" begin
        @test !isnothing(expo_func)
    end

    expected_value = 0.05000045000202501
    tolerance = 1e-3
    @testset "Check expo_func accuracy (low Rt value)" begin
        diff = abs(expo_func - expected_value)
        @test diff <= tolerance
    end

end
