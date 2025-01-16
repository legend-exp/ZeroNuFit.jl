using Pkg
Pkg.activate(".") # activate the environment
Pkg.instantiate() # instantiate the environment
include("../../src/ZeroNuFit.jl")
using .ZeroNuFit
include("../../main.jl")
include("../../src/utils.jl")

@testset "" begin

    @info "Testing function to flag if an event is contained or not in the fit ranges (function 'event_is_contained' in src/utils.jl)"
    fit_range = [[1930.0, 2098.511], [2108.511, 2113.513], [2123.513, 2190.0]]

    # contained event
    event = 1950.0
    flag = nothing
    try
        flag = ZeroNuFit.event_is_contained(event, fit_range)
    catch e
        @error "Error in event_is_contained: $e"
        throw(e)
    end

    @testset "Check flag is valid" begin
        @test !isnothing(flag)
    end
    @test flag == true

    # event outside fit window
    event = 150.0
    flag = ZeroNuFit.event_is_contained(event, fit_range)
    @test flag == false
end
