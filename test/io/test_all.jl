using Test

Test.@testset "likelihood" begin
    include("test_get_partitions.jl")
    include("test_get_events.jl")
    include("test_get_partition_event_index.jl")
    include("test_get_partitions_events.jl")
end
