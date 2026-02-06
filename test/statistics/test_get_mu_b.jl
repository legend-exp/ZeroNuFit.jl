using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("../../src/ZeroNuFit.jl")
using .ZeroNuFit
using Test

@testset "test_get_mu_b" begin

    @info "Testing function to get expected background counts (function 'get_mu_b' in src/likelihood.jl)"

    # Test case 1: Using ckky units (counts/keV/kg/yr)
    @testset "ckky units" begin
        deltaE = 50.0  # keV
        exposure = 100.0  # kg·yr
        bkg_index = 0.01  # counts/keV/kg/yr
        reso = 2.0  # Not used for ckky
        bkg_units = "ckky"

        result = ZeroNuFit.Likelihood.get_mu_b(deltaE, exposure, bkg_index, reso, bkg_units)
        expected = deltaE * exposure * bkg_index  # 50 * 100 * 0.01 = 50

        @test result ≈ expected rtol=1e-10
    end

    # Test case 2: Using cFty units (counts/FWHM/t/yr)
    @testset "cFty units" begin
        deltaE = 50.0  # keV
        exposure = 100.0  # kg·yr
        bkg_index = 1.0  # counts/FWHM/t/yr
        reso = 2.0  # keV (FWHM = 2.0 * 2.355)
        bkg_units = "cFty"

        result = ZeroNuFit.Likelihood.get_mu_b(deltaE, exposure, bkg_index, reso, bkg_units)
        fwhm = reso * 2.355
        expected = deltaE * exposure * (bkg_index / fwhm / 1000)

        @test result ≈ expected rtol=1e-10
    end

    # Test case 3: Zero background index (ckky)
    @testset "Zero background (ckky)" begin
        deltaE = 50.0
        exposure = 100.0
        bkg_index = 0.0
        reso = 2.0
        bkg_units = "ckky"

        result = ZeroNuFit.Likelihood.get_mu_b(deltaE, exposure, bkg_index, reso, bkg_units)
        @test result ≈ 0.0 atol=1e-10
    end

    # Test case 4: Scaling with deltaE
    @testset "Scaling with deltaE (ckky)" begin
        exposure = 100.0
        bkg_index = 0.01
        reso = 2.0
        bkg_units = "ckky"

        result1 = ZeroNuFit.Likelihood.get_mu_b(50.0, exposure, bkg_index, reso, bkg_units)
        result2 = ZeroNuFit.Likelihood.get_mu_b(100.0, exposure, bkg_index, reso, bkg_units)

        # Result should scale linearly with deltaE
        @test result2 ≈ 2.0 * result1 rtol=1e-10
    end

    # Test case 5: Scaling with exposure
    @testset "Scaling with exposure (ckky)" begin
        deltaE = 50.0
        bkg_index = 0.01
        reso = 2.0
        bkg_units = "ckky"

        result1 = ZeroNuFit.Likelihood.get_mu_b(deltaE, 100.0, bkg_index, reso, bkg_units)
        result2 = ZeroNuFit.Likelihood.get_mu_b(deltaE, 200.0, bkg_index, reso, bkg_units)

        # Result should scale linearly with exposure
        @test result2 ≈ 2.0 * result1 rtol=1e-10
    end

    # Test case 6: cFty with different resolution
    @testset "cFty resolution dependence" begin
        deltaE = 50.0
        exposure = 100.0
        bkg_index = 1.0
        bkg_units = "cFty"

        result1 = ZeroNuFit.Likelihood.get_mu_b(deltaE, exposure, bkg_index, 2.0, bkg_units)
        result2 = ZeroNuFit.Likelihood.get_mu_b(deltaE, exposure, bkg_index, 4.0, bkg_units)

        # Result should be inversely proportional to FWHM (resolution * 2.355)
        @test result2 ≈ 0.5 * result1 rtol=1e-10
    end

end
