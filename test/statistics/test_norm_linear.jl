using Pkg
Pkg.activate(".")
Pkg.instantiate()
using Random
include("../../src/ZeroNuFit.jl")
using .ZeroNuFit
include("../../main.jl")
include("../../src/utils.jl")
include("../../src/constants.jl")

@testset "test_norm_linear" begin

    @info "Testing normalised linear function (function 'norm_linear' in src/fitting.jl)"

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
    linear_func = nothing
    try
        linear_func = ZeroNuFit.norm_linear(x, p, b_name, fit_range)
    catch e
        @error "Error in 'norm_linear' evaluation: $e"
        throw(e)
    end

    @testset "Check linear_func is valid" begin
        @test !isnothing(linear_func)
    end

    expected_value = 0.08333333333333333
    tolerance = 1e-3
    @testset "Check linear_func accuracy" begin
        diff = abs(linear_func - expected_value)
        @test diff <= tolerance
    end

end
