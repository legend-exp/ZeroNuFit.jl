using Pkg
Pkg.activate(".")
Pkg.instantiate()
using Random
include("../../src/ZeroNuFit.jl")
using .ZeroNuFit

@testset "test_get_signal_prior_info" begin

    @info "Testing function to retrieve signal prior info (function 'get_signal_prior_info' in src/likelihood.jl)"

    # flat S prior (default)
    config = Dict("signal" => Dict("prior" => "uniform", "upper_bound" => 10))

    # only B fit
    bkg_only = true
    sqrt_prior = nothing
    s_max = nothing
    try
        sqrt_prior, s_max = ZeroNuFit.Likelihood.get_signal_prior_info(bkg_only, config)
    catch e
        @error "Error in 'get_signal_prior_info' evaluation: $e"
        throw(e)
    end

    @test sqrt_prior == false
    @test s_max == nothing

    # B+S fit
    bkg_only = false
    sqrt_prior = nothing
    s_max = nothing
    try
        sqrt_prior, s_max = ZeroNuFit.Likelihood.get_signal_prior_info(bkg_only, config)
    catch e
        @error "Error in 'get_signal_prior_info' evaluation: $e"
        throw(e)
    end

    @test sqrt_prior == false
    @test s_max == nothing

    # 1/sqrt(S) prior
    config = Dict("signal" => Dict("prior" => "sqrt", "upper_bound" => 10))

    bkg_only = false
    sqrt_prior = nothing
    s_max = nothing
    try
        sqrt_prior, s_max = ZeroNuFit.Likelihood.get_signal_prior_info(bkg_only, config)
    catch e
        @error "Error in 'get_signal_prior_info' evaluation: $e"
        throw(e)
    end

    @test sqrt_prior == true
    @test s_max == 10.0

end
