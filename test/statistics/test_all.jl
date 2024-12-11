using Test

Test.@testset "likelihood" begin
    include("test_get_signal_pdf.jl")
    #include("test_get_mu_s_b.jl")
    #include("test_build_likelihood_zero_obs_evts.jl")
    #include("test_build_likelihood_per_partition.jl")
end