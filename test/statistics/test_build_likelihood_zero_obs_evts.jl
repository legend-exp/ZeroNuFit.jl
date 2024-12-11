using Pkg
Pkg.activate(".")
Pkg.instantiate()
using Random
using Test
include("../../src/ZeroNuFit.jl")
using .ZeroNuFit
include("../../main.jl")
include("../../src/utils.jl")

@testset "test_likelihood_one_partition_0_evts" begin
    
    @info "Testing likelihood for one partition but 0 events (function 'build_likelihood_zero_obs_evts' in src/build_likelihood_zero_obs_evts.jl)"
    
    events = [[]]
    
    fit_ranges = Dict("l200a" => [[1930.0, 2098.511], [2108.511, 2113.513], [2123.513, 2190.0]])
    
    # define one partition
    part_event_index = [1]
    partitions = Table(experiment=Array(["L200"]),
                fit_group=Array(["l200a"]),
                bkg_name=Array([:B_l200a_all]),
                signal_name=Array([:gaussian]),
                energy_reso_name=Array([:αr_all]),
                energy_bias_name=Array([:αb_all]),
                eff_name=Array([:αe_all]),
                detector=Array(["V09374A"]),
                part_name=Array(["part0002"]),
                start_ts=Array([1686522477]),
                end_ts=Array([1690190843]),
                eff_tot=Array([0.61102895307898]),
                eff_tot_sigma=Array([0.02786287862866617]),
                width=Array([1.3342225765163704]),
                width_sigma=Array([0.012426367618095815]),
                exposure=Array([0.29072317286485666]),
                bias=Array([-0.04124449186827928]),
                bias_sigma=Array([0.09884528165838884]),
                frac=Array([nothing]),
                tau=Array([nothing]),
                sigma=Array([nothing]))
    
    settings=Dict()
    settings[:energy_scale_fixed] = true
    settings[:energy_scale_correlated] = false
    settings[:eff_fixed] = true
    settings[:eff_correlated] = true
    settings[:bkg_only] = false
    
    sqrt_prior = false
    s_max = nothing
    bkg_shape = :uniform
    p = (S = 100, αe_all = 0.1, ω = [1.1], 𝛥 = [0.1], B_l200a_all = 2E-4)
    
    ll_value = nothing
    try
        ll_value = ZeroNuFit.build_likelihood_zero_obs_evts(partitions[1], p, settings, fit_ranges[partitions[1].fit_group])
    catch e
        @error "Error in 'build_likelihood_zero_obs_evts' evaluation: $e"
        throw(e)
    end

    @testset "Check ll_value is valid" begin
        @test !isnothing(ll_value)
    end
    
    expected_value = -0.11206788736159046
    tolerance = 1e-3      
    @testset "Check ll_value accuracy" begin
        diff = abs(ll_value - expected_value)
        @test diff <= tolerance
    end

end
