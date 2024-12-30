using Pkg
Pkg.activate(".")
Pkg.instantiate()
using Random
include("../../src/ZeroNuFit.jl")
using .ZeroNuFit
include("../../main.jl")

@testset "test_get_bkg_info" begin
    
    @info "Testing gaussian signal for one partition and one close event @ 2039.05 keV (function 'get_bkg_info' in src/fitting.jl)"
    
    # default modeling, i.e. flat
    config = Dict(
        "bkg" => Dict()
    )
    
    bkg_shape = nothing
    bkg_shape_pars = nothing
    try
        bkg_shape,bkg_shape_pars = ZeroNuFit.get_bkg_info(config)
    catch e
        @error "Error in 'get_bkg_info' evaluation: $e"
        throw(e)
    end

    expected_value = :uniform
    @testset "Check bkg_shape accuracy" begin
        @test bkg_shape == expected_value
    end
    expected_value = nothing
    @testset "Check bkg_shape_pars accuracy" begin
        @test bkg_shape_pars == expected_value
    end
    
    # linear modeling
    config = Dict("bkg" => Dict("shape" => Dict("name" => "linear", "pars" => Dict("slope" => [-10, 10]))))
    
    bkg_shape = nothing
    bkg_shape_pars = nothing
    try
        bkg_shape,bkg_shape_pars = ZeroNuFit.get_bkg_info(config)
    catch e
        @error "Error in 'get_bkg_info' evaluation: $e"
        throw(e)
    end
    
    expected_value = :linear
    @testset "Check bkg_shape accuracy" begin
        @test bkg_shape == expected_value
    end
    expected_value = Dict("slope" => [-10, 10])
    @testset "Check bkg_shape_pars accuracy" begin
        @test bkg_shape_pars == expected_value
    end
    
    # exponential modeling
    config = Dict("bkg" => Dict("shape" => Dict("name" => "exponential", "pars" => Dict("slope" => [-5, 5]))))
    
    bkg_shape = nothing
    bkg_shape_pars = nothing
    try
        bkg_shape,bkg_shape_pars = ZeroNuFit.get_bkg_info(config)
    catch e
        @error "Error in 'get_bkg_info' evaluation: $e"
        throw(e)
    end
    
    expected_value = :exponential
    @testset "Check bkg_shape accuracy" begin
        @test bkg_shape == expected_value
    end
    expected_value = Dict("slope" => [-5, 5])
    @testset "Check bkg_shape_pars accuracy" begin
        @test bkg_shape_pars == expected_value
    end
    
end
