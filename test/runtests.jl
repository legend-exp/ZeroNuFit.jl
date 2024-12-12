using Test

@info "Running tests with $(Base.Threads.nthreads()) Julia threads active."

import Logging
Logging.global_logger(Logging.ConsoleLogger(stdout, Logging.Info))

Test.@testset "Package ZeroNuFit" begin
    include("io/test_all.jl")
    include("statistics/test_all.jl")
end