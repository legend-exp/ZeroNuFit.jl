using Pkg
Pkg.activate(".") # activate the environment
Pkg.instantiate() # instantiate the environment
include("../../src/ZeroNuFit.jl")
using .ZeroNuFit
include("../../main.jl")

@testset "test_get_range" begin
    
    @info "Testing function to retrieve lower and upper range edges (function 'get_range' in src/fitting.jl)"
    
    # 1 entry in fit ranges
    fit_ranges = [[1930.0, 1950.0]]
    range_l = nothing 
    range_h = nothing
    try
        range_l, range_h = ZeroNuFit.get_range(fit_ranges)
    catch e
        @error "Error in get_range: $e"
        throw(e)
    end
    
    @testset "Check range_l is valid" begin
        @test !isnothing(range_l)
    end
    @testset "Check range_h is valid" begin
        @test !isnothing(range_h)
    end
    
    expected_events=[1930.0]
    @testset "Check events accuracy" begin
        @test range_l == expected_events
    end
    expected_events=[1950.0]
    @testset "Check events accuracy" begin
        @test range_h == expected_events
    end
    
    # more entries in fit ranges
    fit_ranges = [[1930.0, 1950.0], [2100.0,2200.0], [1450.0, 1500.0]]
    range_l = nothing 
    range_h = nothing
    try
        range_l, range_h = ZeroNuFit.get_range(fit_ranges)
    catch e
        @error "Error in get_range: $e"
        throw(e)
    end
    
    @testset "Check range_l is valid" begin
        @test !isnothing(range_l)
    end
    @testset "Check range_h is valid" begin
        @test !isnothing(range_h)
    end
    
    expected_events=[1450.0, 1930.0, 2100.0]
    @testset "Check events accuracy" begin
        @test range_l == expected_events
    end
    expected_events=[1500.0, 1950.0, 2200.0]
    @testset "Check events accuracy" begin
        @test range_h == expected_events
    end
end
    
