using Pkg
Pkg.activate(".") # activate the environment
Pkg.instantiate() # instantiate the environment
include("../../src/ZeroNuFit.jl")
using .ZeroNuFit

struct Sample
    v::Dict{Symbol,Any}
    weight::Int
end

@testset "test_get_par_posterior " begin

    @info "Testing function to retrieve samples of a posterior pdf (function 'get_par_posterior ' in src/utils.jl)"

    # testing when idx is nothing
    @testset "Test get_par_posterior function" begin
        sample1 = Sample(Dict(:alpha => [1, 2, 3], :beta => [4, 5]), 2)
        sample2 = Sample(Dict(:alpha => [6, 7, 8], :beta => [9, 10]), 1)
        samples = [sample1, sample2]

        # test when idx is nothing (should return the entire value of par)
        result = ZeroNuFit.Utils.get_par_posterior(samples, :alpha)
        @test result == [1, 2, 3, 1, 2, 3, 6, 7, 8]
    end

    #testing when idx is a specific index (should return the indexed value)
    @testset "Test get_par_posterior function with idx" begin
        sample1 = Sample(Dict(:alpha => [1, 2, 3], :beta => [4, 5]), 2)
        sample2 = Sample(Dict(:alpha => [6, 7, 8], :beta => [9, 10]), 1)
        samples = [sample1, sample2]

        # test when idx is 2 (should return the 2nd element of :alpha)
        result = ZeroNuFit.Utils.get_par_posterior(samples, :alpha, idx = 2)
        @test result == [2, 2, 7]
    end

    # empty sample list
    @testset "Edge case with empty samples" begin
        samples = []

        # test with empty sample list (should return an empty result)
        result = ZeroNuFit.Utils.get_par_posterior(samples, :alpha)
        @test result == []
    end

    # case where a parameter doesn't exist in the sample
    @testset "Test missing parameter" begin
        sample1 = Sample(Dict(:alpha => [1, 2, 3], :beta => [4, 5]), 1)
        samples = [sample1]

        # test when the parameter doesn't exist (should raise an error or return an empty list)
        try
            result = ZeroNuFit.Utils.get_par_posterior(samples, :gamma)
            @test result == []
        catch e
            @test true
        end
    end

end
