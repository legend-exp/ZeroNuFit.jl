using Pkg
Pkg.activate(".")
Pkg.instantiate()
using Random
include("../../src/ZeroNuFit.jl")
using .ZeroNuFit
include("../../main.jl")
include("../../src/utils.jl")
include("../../src/constants.jl")

Base.exit(code::Int) = throw(ArgumentError("exit code $code"))

@testset "test_get_signal_pdf" begin
    
    @info "Testing gaussian signal for one partition and one close event @ 2039.05 keV (function 'get_signal_pdf' in src/likelihood.jl)"
    
    event = 2039.05
    
    # gaussian pdf
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
                width=Array([1.3]),
                width_sigma=Array([0.01]),
                exposure=Array([1]),
                bias=Array([0.1]),
                bias_sigma=Array([0.01]),
                frac=Array([nothing]),
                tau=Array([nothing]),
                sigma=Array([nothing]))
    
    signal_pdf = nothing
    try
        signal_pdf = ZeroNuFit.get_signal_pdf(event,constants.Qbb,partitions[1])
    catch e
        @error "Error in 'get_signal_pdf' evaluation: $e"
        throw(e)
    end

    @testset "Check signal_pdf is valid" begin
        @test !isnothing(signal_pdf)
    end
    
    expected_value = 0.3061441384108178
    tolerance = 1e-5      
    @testset "Check signal_pdf accuracy" begin
        diff = abs(signal_pdf - expected_value)
        @test diff <= tolerance
    end
    
    
    # gaussian + low-E tail signal pdf
    partitions = Table(experiment=Array(["L200"]),
                fit_group=Array(["l200a"]),
                bkg_name=Array([:B_l200a_all]),
                signal_name=Array([:gaussian_plus_lowEtail]),
                energy_reso_name=Array([:αr_all]),
                energy_bias_name=Array([:αb_all]),
                eff_name=Array([:αe_all]),
                detector=Array(["V09374A"]),
                part_name=Array(["part0002"]),
                start_ts=Array([1686522477]),
                end_ts=Array([1690190843]),
                eff_tot=Array([0.5]),
                eff_tot_sigma=Array([0.01]),
                width=Array([1.0]),
                width_sigma=Array([0.05]),
                exposure=Array([1]),
                bias=Array([0.1]),
                bias_sigma=Array([0.01]),
                frac=Array([0.1]),
                tau=Array([1.5]),
                sigma=Array([1.1]))
    
    
    signal_pdf = nothing
    try
        signal_pdf = ZeroNuFit.get_signal_pdf(event,constants.Qbb,partitions[1])
    catch e
        @error "Error in 'get_signal_pdf' evaluation: $e"
        throw(e)
    end

    @testset "Check signal_pdf is valid" begin
        @test !isnothing(signal_pdf)
    end
    
    expected_value = 0.3445363167418973
    tolerance = 1e-5      
    @testset "Check signal_pdf accuracy" begin
        diff = abs(signal_pdf - expected_value)
        @test diff <= tolerance
    end
    
    
    # not-existing signal shape
    partitions = Table(experiment=Array(["L200"]),
                fit_group=Array(["l200a"]),
                bkg_name=Array([:B_l200a_all]),
                signal_name=Array([:fancy_gaussian]),
                energy_reso_name=Array([:αr_all]),
                energy_bias_name=Array([:αb_all]),
                eff_name=Array([:αe_all]),
                detector=Array(["V09374A"]),
                part_name=Array(["part0002"]),
                start_ts=Array([1686522477]),
                end_ts=Array([1690190843]),
                eff_tot=Array([0.5]),
                eff_tot_sigma=Array([0.01]),
                width=Array([1.0]),
                width_sigma=Array([0.05]),
                exposure=Array([1]),
                bias=Array([0.1]),
                bias_sigma=Array([0.01]),
                frac=Array([0.1]),
                tau=Array([1.5]),
                sigma=Array([1.1]))
    @test_throws ArgumentError ZeroNuFit.get_signal_pdf(event, constants.Qbb, partitions[1])
end
