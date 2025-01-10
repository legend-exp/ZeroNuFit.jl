using Pkg
Pkg.activate(".") # activate the environment
Pkg.instantiate() # instantiate the environment
include("../../src/ZeroNuFit.jl")
using .ZeroNuFit
include("../../main.jl")
include("../../src/utils.jl")

@testset "test_get_energy_scale_pars" begin
    
    @info "Testing function to retrieve resolution and bias (function 'get_energy_scale_pars' in src/utils.jl)"
    
    # fixed energy scale OR no events in the partition
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
                eff_tot=Array([0.5]),
                eff_tot_sigma=Array([0.01]),
                width=Array([0.5]),
                width_sigma=Array([0.01]),
                exposure=Array([1]),
                bias=Array([0.5]),
                bias_sigma=Array([0.01]),
                frac=Array([nothing]),
                tau=Array([nothing]),
                sigma=Array([nothing]))
    part_event_index = [1]
    p = (S = 100, αe_all = 0.5, αb_all=0.5, αr_all=0.5, ω = [1.0], 𝛥 = [1.0], B_l200a_all = 2E-4)
    settings =Dict(
        :energy_scale_fixed => true,
        :energy_scale_correlated => false
    )
    reso = nothing
    bias = nothing
    # fixed nuisance, 1 event in the partition
    try
        reso,bias = ZeroNuFit.get_energy_scale_pars(partitions[1],p,settings,1)
    catch e
        @error "Error in get_energy_scale_pars: $e"
        throw(e)
    end
    
    @testset "Check resolution/bias is valid" begin
        @test !isnothing(reso)
        @test !isnothing(bias)
    end
    
    expected_reso = 0.5
    expected_bias = 0.5
    @testset "Check resolution/bias accuracy [fixed nuisance, 1 event in the partition]" begin
        @test reso == expected_reso
        @test bias == expected_bias
    end
    
    # fixed nuisance, 0 event in the partition
    reso = nothing
    bias = nothing
    reso,bias = ZeroNuFit.get_energy_scale_pars(partitions[1],p,settings,0)
    expected_reso = 0.5
    expected_bias = 0.5
    @testset "Check resolution/bias accuracy [fixed nuisance, 0 event in the partition]" begin
        @test reso == expected_reso
        @test bias == expected_bias
    end
    
    # not fixed nuisance, 0 event in the partition
    reso = nothing
    bias = nothing
    settings =Dict(
        :energy_scale_fixed => false,
        :energy_scale_correlated => false
    )
    reso,bias = ZeroNuFit.get_energy_scale_pars(partitions[1],p,settings,0)
    expected_reso = 0.5
    expected_bias = 0.5
    @testset "Check resolution/bias accuracy [not fixed nuisance, 0 event in the partition]" begin
        @test reso == expected_reso
        @test bias == expected_bias
    end
    
    # correlated energy scale, 1 event in the partition
    reso = nothing
    bias = nothing
    settings =Dict(
        :energy_scale_fixed => false,
        :energy_scale_correlated => true
    )
    reso,bias = ZeroNuFit.get_energy_scale_pars(partitions[1],p,settings,1)
    energy_reso_group = partitions[1].energy_reso_name
    energy_bias_group = partitions[1].energy_bias_name
    expected_reso = partitions[1].width+p[energy_reso_group]*partitions[1].width_sigma
    expected_bias = partitions[1].bias+p[energy_bias_group]*partitions[1].bias_sigma
    @testset "Check resolution/bias accuracy [correlated energy scale, 1 event in the partition]" begin
        @test reso == expected_reso
        @test bias == expected_bias
    end
    
    # correlated energy scale, 0 event in the partition
    reso = nothing
    bias = nothing
    settings =Dict(
        :energy_scale_fixed => false,
        :energy_scale_correlated => true
    )
    reso,bias = ZeroNuFit.get_energy_scale_pars(partitions[1],p,settings,0)
    expected_reso = 0.5
    expected_bias = 0.5
    @testset "Check resolution/bias accuracy [correlated energy scale, 0 event in the partition]" begin
        @test reso == expected_reso
        @test bias == expected_bias
    end
    
    # uncorrelated energy scale, 1 event in the partition
    reso = nothing
    bias = nothing
    settings =Dict(
        :energy_scale_fixed => false,
        :energy_scale_correlated => false
    )
    reso,bias = ZeroNuFit.get_energy_scale_pars(partitions[1],p,settings,1)
    expected_reso = 1.0
    expected_bias = 1.0
    @testset "Check resolution/bias accuracy [uncorrelated energy scale, 1 event in the partition]" begin
        @test reso == expected_reso
        @test bias == expected_bias
    end
    
    # uncorrelated energy scale, 0 event in the partition
    reso = nothing
    bias = nothing
    settings =Dict(
        :energy_scale_fixed => false,
        :energy_scale_correlated => false
    )
    reso,bias = ZeroNuFit.get_energy_scale_pars(partitions[1],p,settings,0)
    expected_reso = 0.5
    expected_bias = 0.5
    @testset "Check resolution/bias accuracy [uncorrelated energy scale, 0 event in the partition]" begin
        @test reso == expected_reso
        @test bias == expected_bias
    end
    
end
    
