using Test

Test.@testset "io" begin
    include("test_get_settings.jl")
    include("test_get_partitions.jl")
    include("test_get_events.jl")
    include("test_get_partition_event_index.jl")
    include("test_get_partitions_events.jl")
    include("test_get_corr_info.jl")
    include("test_get_signal_prior_info.jl")
    include("test_get_energy_scale_pars.jl")
    include("test_get_efficiency.jl")
end
