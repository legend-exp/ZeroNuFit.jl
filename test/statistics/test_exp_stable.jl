using Pkg
Pkg.activate(".")
Pkg.instantiate()
using Random
include("../../src/ZeroNuFit.jl")
using .ZeroNuFit

@testset "test_exp_stable" begin

    @info "Testing exponential function for low/high x (function 'exp_stable' in src/likelihood.jl)"

    # small values where Taylor expansion is valid
    @test abs(ZeroNuFit.Likelihood.exp_stable(0.0) - 1.0) < 1E-12
    @test abs(ZeroNuFit.Likelihood.exp_stable(1E-7) - (1 + 1E-7 + (1E-7)^2 / 2 + (1E-7)^3 / 6)) < 1E-12
    @test abs(ZeroNuFit.Likelihood.exp_stable(-1E-7) - (1 - 1E-7 + (-1E-7)^2 / 2 + (-1E-7)^3 / 6)) <
          1E-12

    # large values where exp(x) is valid
    @test abs(ZeroNuFit.Likelihood.exp_stable(1.0) - exp(1.0)) < 1E-12
    @test abs(ZeroNuFit.Likelihood.exp_stable(-1.0) - exp(-1.0)) < 1E-12

    # zero case
    @test ZeroNuFit.Likelihood.exp_stable(0.0) == 1.0

    # large positive/negative values
    @test abs(ZeroNuFit.Likelihood.exp_stable(10.0) - exp(10.0)) < 1E-12
    @test abs(ZeroNuFit.Likelihood.exp_stable(-10.0) - exp(-10.0)) < 1E-12

    # abs(x) = 1E-6 (transition point between Taylor expansion and exp(x))
    @test abs(ZeroNuFit.Likelihood.exp_stable(1E-6) - exp(1E-6)) < 1E-12
    @test abs(ZeroNuFit.Likelihood.exp_stable(-1E-6) - exp(-1E-6)) < 1E-12


end
