using Pkg
Pkg.activate(".")
Pkg.instantiate()
using Random
include("../../src/ZeroNuFit.jl")
using .ZeroNuFit
include("../../main.jl")
include("../../src/mbb.jl")
include("../../src/constants.jl")

@testset "test_mbb" begin

    @info "Testing function for evaluating mbb (function 'mbb' in src/mbb.jl)"

    # frequentist half-life limit taken from GSTR-20-006
    gerda_T12 = 1.83
    gerda_low_limit = ZeroNuFit.mbb(1 / gerda_T12 * 10, G = constants.phase_space, M = constants.nme_gerda_up)
    gerda_upp_limit = ZeroNuFit.mbb(1 / gerda_T12 * 10, G = constants.phase_space, M = constants.nme_gerda_low)

    expected_prl_low_limit = 79
    expected_prl_upp_limit = 180
    @test expected_prl_low_limit == round(gerda_low_limit)
    @test expected_prl_upp_limit == round(gerda_upp_limit)

end
