using Pkg
Pkg.activate(".")
Pkg.instantiate()
using Random
include("../../src/ZeroNuFit.jl")
using .ZeroNuFit

Base.exit(code::Int) = throw(ArgumentError("exit code $code"))

@testset "test_get_bkg_pdf" begin

    @info "Testing background pdf retrieval (function 'get_bkg_pdf' in src/likelihood.jl)"

    # exponential bkg
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
    bkg = nothing
    try
        bkg = ZeroNuFit.Likelihood.get_bkg_pdf(:exponential, x, p, b_name, fit_range)
    catch e
        @error "Error in 'get_bkg_pdf' evaluation: $e"
        throw(e)
    end

    @testset "Check bkg is valid" begin
        @test !isnothing(bkg)
    end

    expected_value = 0.08506327727340002
    tolerance = 1e-3
    @testset "Check bkg accuracy" begin
        diff = abs(bkg - expected_value)
        @test diff <= tolerance
    end

    # linear bkg
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
    bkg = nothing
    try
        bkg = ZeroNuFit.Likelihood.get_bkg_pdf(:linear, x, p, b_name, fit_range)
    catch e
        @error "Error in 'get_bkg_pdf' evaluation: $e"
        throw(e)
    end

    @testset "Check bkg is valid" begin
        @test !isnothing(bkg)
    end

    expected_value = 0.08333333333333333
    tolerance = 1e-3
    @testset "Check bkg accuracy" begin
        diff = abs(bkg - expected_value)
        @test diff <= tolerance
    end

    # not-existing bkg shape
    @test_throws ArgumentError ZeroNuFit.Likelihood.get_bkg_pdf(
        :fancy_bkg,
        x,
        p,
        b_name,
        fit_range,
    )
end
