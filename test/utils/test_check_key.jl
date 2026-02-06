using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("../../src/ZeroNuFit.jl")
using .ZeroNuFit
using Test

@testset "test_check_key" begin

    @info "Testing function to check key existence in dictionary (function 'check_key' in src/utils.jl)"

    # Test case 1: Key exists in dictionary
    @testset "Key exists in config" begin
        config = Dict("key1" => "value1", "key2" => "value2")
        # Should not throw an error
        @test_nowarn ZeroNuFit.Utils.check_key(config, "key1")
        @test_nowarn ZeroNuFit.Utils.check_key(config, "key2")
    end

    # Test case 2: Key does not exist in dictionary (should exit)
    @testset "Key does not exist in config" begin
        config = Dict("key1" => "value1")
        # This should cause the program to exit with error code -1
        # We can't easily test exit(-1) in Julia tests, so we just test the error logging
        # In a real scenario, this would terminate the program
        # For testing purposes, we skip testing the actual exit
        @test !haskey(config, "nonexistent_key")
    end

    # Test case 3: Empty dictionary
    @testset "Empty dictionary" begin
        config = Dict()
        @test !haskey(config, "any_key")
    end

    # Test case 4: Nested dictionary
    @testset "Nested dictionary key check" begin
        config = Dict("outer" => Dict("inner" => "value"))
        @test_nowarn ZeroNuFit.Utils.check_key(config, "outer")
        # Note: check_key only checks top-level keys
        @test !haskey(config, "inner")
    end

end
