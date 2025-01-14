using Pkg
Pkg.activate(".") # activate the environment
Pkg.instantiate() # instantiate the environment
include("../../src/ZeroNuFit.jl")
using .ZeroNuFit
include("../../main.jl")
include("../../src/utils.jl")

@testset "test_get_efficiency" begin

    @info "Testing function to retrieve efficiency (function 'get_efficiency' in src/utils.jl)"

    partitions = Table(
        experiment = Array(["L200"]),
        fit_group = Array(["l200a"]),
        bkg_name = Array([:B_l200a_all]),
        signal_name = Array([:gaussian]),
        energy_reso_name = Array([:αr_all]),
        energy_bias_name = Array([:αb_all]),
        eff_name = Array([:αe_all]),
        detector = Array(["V09374A"]),
        part_name = Array(["part0002"]),
        start_ts = Array([1686522477]),
        end_ts = Array([1690190843]),
        eff_tot = Array([0.5]),
        eff_tot_sigma = Array([0.01]),
        width = Array([1.3]),
        width_sigma = Array([0.01]),
        exposure = Array([1]),
        bias = Array([0.1]),
        bias_sigma = Array([0.01]),
        frac = Array([nothing]),
        tau = Array([nothing]),
        sigma = Array([nothing]),
    )
    p = (
        S = 100,
        αe_all = 0.1,
        αb_all = 0.5,
        αr_all = 0.5,
        ω = [1.1],
        𝛥 = [0.5],
        B_l200a_all = 2E-4,
        ε = 1,
    )
    settings = Dict(:eff_fixed => true, :eff_correlated => false, :bkg_only => false)
    eff = nothing

    # fixed efficiency, 1 event in the partition
    eff = nothing
    try
        eff = ZeroNuFit.get_efficiency(p, partitions[1], 1, settings)
    catch e
        @error "Error in get_efficiency: $e"
        throw(e)
    end

    @testset "Check eff is valid" begin
        @test !isnothing(eff)
    end

    @testset "Check eff accuracy [fixed efficiency, 1 event in the partition]" begin
        @test eff == 0.5
    end

    # fixed efficiency, 0 event in the partition
    eff = nothing
    eff = ZeroNuFit.get_efficiency(p, partitions[1], 0, settings)
    @testset "Check eff accuracy [fixed efficiency, 0 event in the partition]" begin
        @test eff == 0.5
    end

    # correlated efficiency, 1 event in the partition
    eff = nothing
    settings = Dict(:eff_fixed => false, :eff_correlated => true, :bkg_only => false)
    eff = ZeroNuFit.get_efficiency(p, partitions[1], 1, settings)
    energy_eff_group = partitions[1].eff_name
    expected_eff = partitions[1].eff_tot + p[energy_eff_group] * partitions[1].eff_tot_sigma
    @testset "Check eff accuracy [correlated efficiency, 1 event in the partition]" begin
        @test eff == expected_eff
    end

    # correlated efficiency, 0 event in the partition
    eff = nothing
    settings = Dict(:eff_fixed => false, :eff_correlated => true, :bkg_only => false)
    eff = ZeroNuFit.get_efficiency(p, partitions[1], 0, settings)
    energy_eff_group = partitions[1].eff_name
    expected_eff = partitions[1].eff_tot + p[energy_eff_group] * partitions[1].eff_tot_sigma
    @testset "Check eff accuracy [correlated efficiency, 0 event in the partition]" begin
        @test eff == expected_eff
    end

    # uncorrelated efficiency, 1 event in the partition
    eff = nothing
    settings = Dict(:eff_fixed => false, :eff_correlated => false, :bkg_only => false)
    eff = ZeroNuFit.get_efficiency(p, partitions[1], 1, settings)
    @testset "Check eff accuracy [uncorrelated efficiency, 1 event in the partition]" begin
        @test eff == 1
    end

    # uncorrelated efficiency, 0 event in the partition 
    eff = nothing
    settings = Dict(:eff_fixed => false, :eff_correlated => false, :bkg_only => false)
    eff = ZeroNuFit.get_efficiency(p, partitions[1], 0, settings)
    @testset "Check eff accuracy [uncorrelated efficiency, 0 event in the partition]" begin
        @test eff == 0.5
    end

    # only bkg
    eff = nothing
    settings = Dict(:eff_fixed => false, :eff_correlated => false, :bkg_only => true)
    eff = ZeroNuFit.get_efficiency(p, partitions[1], 0, settings)
    @testset "Check eff accuracy [only bkg]" begin
        @test eff == 0
    end

end
