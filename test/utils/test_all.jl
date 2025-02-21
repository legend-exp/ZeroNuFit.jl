using Test

Test.@testset "likelihood" begin
    include("test_get_range.jl")
    include("test_get_par_posterior.jl")
end
