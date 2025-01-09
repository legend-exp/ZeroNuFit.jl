using Test

Test.@testset "statistics" begin
    include("test_get_signal_pdf.jl")
    include("test_get_bkg_pdf.jl")
    include("test_get_mu_s_b.jl")
    include("test_build_likelihood_zero_obs_evts.jl")
    include("test_build_likelihood_per_partition.jl")
    include("test_inverse_uniform_cdf.jl")
    include("test_generate_disjoint_uniform_samples.jl")
    include("test_get_bkg_info.jl")
    include("test_norm_linear.jl")
    include("test_exp_stable.jl")
    include("test_norm_exponential.jl")
end