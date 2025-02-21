### runtests.jl
#
# Authors: Sofia Calgaro, Toby Dixon
# 
###
using Test

@info "Running tests with $(Base.Threads.nthreads()) Julia threads active."

import Logging
Logging.global_logger(Logging.ConsoleLogger(stdout, Logging.Info))

Test.@testset "Package ZeroNuFit" begin
    include("io/test_all.jl")
    include("plotting/test_all.jl")
    include("statistics/test_all.jl")
    include("utils/test_all.jl")
end
