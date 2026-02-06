using Test

Test.@testset "utils" begin
    include("test_get_range.jl")
    include("test_get_par_posterior.jl")
    include("test_check_key.jl")
    include("test_get_deltaE.jl")
    include("test_get_corr_info.jl")
end
