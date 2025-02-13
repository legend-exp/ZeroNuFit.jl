using Pkg
Pkg.activate(".")
Pkg.instantiate()
using Random
include("../../src/ZeroNuFit.jl")
using .ZeroNuFit

@testset "test_generate_disjoint_uniform_samples" begin

    @info "Testing function for getting random 'n' energies (function 'generate_disjoint_uniform_samples' in src/utils.jl)"

    # fixed seed (1 energy, 1 range)
    one_energy = nothing
    fit_range = [[1920.0, 1930.0]]
    try
        one_energy = ZeroNuFit.Utils.generate_disjoint_uniform_samples(1, fit_range; seed = 123)
    catch e
        @error "Error in 'generate_disjoint_uniform_samples' evaluation: $e"
        throw(e)
    end

    @testset "Check one_energy is valid (fixed seed)" begin
        @test !isnothing(one_energy)
    end

    expected_value = [1925.212137955354]
    @testset "Check energy accuracy (fixed seed)" begin
        @test one_energy == expected_value
    end

    function check_in_ranges(energies, fit_range)
        return [any(e >= r[1] && e <= r[2] for r in fit_range) for e in energies]
    end

    @testset "Check energy containment in fit range (fixed seed)" begin
        contained = check_in_ranges(one_energy, fit_range)
        @test all(contained)
    end

    # random seed generator (1 energy, 1 range)
    one_energy = nothing
    try
        one_energy = ZeroNuFit.Utils.generate_disjoint_uniform_samples(1, fit_range)
    catch e
        @error "Error in 'generate_disjoint_uniform_samples' evaluation: $e"
        throw(e)
    end

    @testset "Check one_energy is valid (random seed)" begin
        @test !isnothing(one_energy)
    end

    expected_value = [1925.212137955354]
    @testset "Check energy accuracy (random seed, 1 energy)" begin
        @test one_energy != expected_value
    end

    @testset "Check energy containment in fit range (random seed, 1 energy)" begin
        contained = check_in_ranges(one_energy, fit_range)
        @test all(contained)
    end

    # random seed generator (more energies, 1 range)
    more_energies = nothing
    try
        more_energies = ZeroNuFit.Utils.generate_disjoint_uniform_samples(10, fit_range)
    catch e
        @error "Error in 'generate_disjoint_uniform_samples' evaluation: $e"
        throw(e)
    end

    @testset "Check more_energies is valid (random seed, more energies)" begin
        @test !isnothing(more_energies)
    end

    @testset "Check energy containment in fit range (random seed, more energies)" begin
        contained = check_in_ranges(more_energies, fit_range)
        @test all(contained)
    end

    # random seed generator (more energies, more ranges)
    more_energies = nothing
    fit_range = [[1920.0, 1930.0], [1970.0, 1980.0], [2100.0, 2250.0]]
    try
        more_energies = ZeroNuFit.Utils.generate_disjoint_uniform_samples(10, fit_range)
    catch e
        @error "Error in 'generate_disjoint_uniform_samples' evaluation: $e"
        throw(e)
    end

    @testset "Check more_energies is valid (random seed, more energies, more ranges)" begin
        @test !isnothing(more_energies)
    end

    @testset "Check energy containment in fit range (random seed, more energies, more ranges)" begin
        contained = check_in_ranges(more_energies, fit_range)
        @test all(contained)
    end

end
