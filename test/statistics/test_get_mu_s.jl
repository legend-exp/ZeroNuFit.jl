using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("../../src/ZeroNuFit.jl")
using .ZeroNuFit
using Test

@testset "test_get_mu_s" begin

    @info "Testing function to get expected signal counts (function 'get_mu_s' in src/likelihood.jl)"

    # Test case 1: Basic calculation
    @testset "Basic signal count calculation" begin
        exposure = 100.0  # kg·yr
        eff = 0.8  # 80% efficiency
        signal = 1.0  # Signal strength in units of 1e-27
        
        result = ZeroNuFit.Likelihood.get_mu_s(exposure, eff, signal)
        
        # Expected: log(2) * N_A * exposure * eff * signal * sig_units / m_76
        N_A = ZeroNuFit.Constants.N_A
        m_76 = ZeroNuFit.Constants.m_76
        sig_units = ZeroNuFit.Constants.sig_units
        expected = log(2) * N_A * exposure * eff * signal * sig_units / m_76
        
        @test result ≈ expected rtol=1e-10
    end

    # Test case 2: Zero efficiency
    @testset "Zero efficiency" begin
        exposure = 100.0
        eff = 0.0
        signal = 1.0
        
        result = ZeroNuFit.Likelihood.get_mu_s(exposure, eff, signal)
        @test result ≈ 0.0 atol=1e-10
    end

    # Test case 3: Zero signal
    @testset "Zero signal" begin
        exposure = 100.0
        eff = 0.8
        signal = 0.0
        
        result = ZeroNuFit.Likelihood.get_mu_s(exposure, eff, signal)
        @test result ≈ 0.0 atol=1e-10
    end

    # Test case 4: Different exposure values
    @testset "Scaling with exposure" begin
        eff = 0.8
        signal = 1.0
        
        result1 = ZeroNuFit.Likelihood.get_mu_s(100.0, eff, signal)
        result2 = ZeroNuFit.Likelihood.get_mu_s(200.0, eff, signal)
        
        # Result should scale linearly with exposure
        @test result2 ≈ 2.0 * result1 rtol=1e-10
    end

    # Test case 5: Different efficiency values
    @testset "Scaling with efficiency" begin
        exposure = 100.0
        signal = 1.0
        
        result1 = ZeroNuFit.Likelihood.get_mu_s(exposure, 0.5, signal)
        result2 = ZeroNuFit.Likelihood.get_mu_s(exposure, 1.0, signal)
        
        # Result should scale linearly with efficiency
        @test result2 ≈ 2.0 * result1 rtol=1e-10
    end

end
