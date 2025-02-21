using Pkg
Pkg.activate(".")
Pkg.instantiate()
using Random
include("../../src/ZeroNuFit.jl")
using .ZeroNuFit
using IntervalSets
using Distributions
using JSON
using LaTeXStrings

@testset "test_build_prior" begin

    @info "Testing retrieval of signal & background pdfs (function 'build_prior' in src/likelihood.jl)"
    present_dir = @__DIR__

    config = Dict(
        "bkg_only" => false,
        "bkg" => Dict(
            "correlated" => Dict("range" => "none", "mode" => "none"),
            "prior" => "uniform",
            "upper_bound" => 0.1,
        ),
        "signal" => Dict("prior" => "uniform", "upper_bound" => 1000),
        "events" => [joinpath(present_dir, "../inputs/events_test.json")],
        "partitions" => [joinpath(present_dir, "../inputs/partitions_test.json")],
        "output_path" => "tests",
        "bat_fit" => Dict("nsteps" => 10000.0, "nchains" => 4),
        "nuisance" => Dict(
            "efficiency" => Dict("fixed" => true, "correlated" => true),
            "energy_bias" => Dict("fixed" => true, "correlated" => false),
            "energy_res" => Dict("fixed" => true, "correlated" => false),
        ),
        "plot" => Dict(
            "bandfit_and_data" => false,
            "alpha" => 0.3,
            "fit_and_data" => false,
            "scheme" => "blue",
        ),
    )

    # fixed efficiency & energy scale parameters
    partitions, fit_ranges = ZeroNuFit.Utils.get_partitions(config)
    events = ZeroNuFit.Utils.get_events(config["events"][1], partitions)
    part_event_index = ZeroNuFit.Utils.get_partition_event_index(events, partitions)
    settings = ZeroNuFit.Utils.get_settings(config)
    bkg_name =
        JSON.parsefile(config["partitions"][1])["fit_groups"]["all_phase_II"]["bkg_name"]
    priors = nothing
    pretty_names = nothing
    nuisance_info = nothing
    try
        priors, pretty_names, nuisance_info =
            ZeroNuFit.Likelihood.build_prior(partitions, part_event_index, config, settings)
    catch e
        @error "Error in build_prior: $e"
        throw(e)
    end
    # S & B priors
    @test Uniform(0, 1000) == priors.S
    @test Uniform(0, 0.1) == priors.B_gerda_all_pII
    # pretty names (just check it once, these are unique)
    @test string("S [") * L"10^{-27}" * string("yr") * L"^{-1}" * string("]") ==
          pretty_names[:S]
    @test string(bkg_name) * " [cts/keV/kg/yr]" == pretty_names[Symbol(bkg_name)]
    @test L"\alpha_{\varepsilon}" == pretty_names[:α]
    @test L"\alpha_{r}" == pretty_names[:αr]
    @test L"\alpha_{b}" == pretty_names[:αb]
    # nuisance_info
    @test [] == nuisance_info["α"]
    @test [] == nuisance_info["αr"]
    @test [] == nuisance_info["αb"]
    @test [] == nuisance_info["γ"]
    @test [] == nuisance_info["ε"]
    @test [] == nuisance_info["ω"]
    @test [] == nuisance_info["𝛥"]

    # varying efficiency & energy scale parameters (uncorrelated)
    config["nuisance"]["efficiency"]["fixed"] = false
    config["nuisance"]["energy_bias"]["fixed"] = false
    config["nuisance"]["energy_res"]["fixed"] = false
    partitions, fit_ranges = ZeroNuFit.Utils.get_partitions(config)
    events = ZeroNuFit.Utils.get_events(config["events"][1], partitions)
    part_event_index = ZeroNuFit.Utils.get_partition_event_index(events, partitions)
    settings = ZeroNuFit.Utils.get_settings(config)
    bkg_name =
        JSON.parsefile(config["partitions"][1])["fit_groups"]["all_phase_II"]["bkg_name"]
    priors = nothing
    pretty_names = nothing
    nuisance_info = nothing
    priors, pretty_names, nuisance_info =
        ZeroNuFit.Likelihood.build_prior(partitions, part_event_index, config, settings)

    # priors
    @test Uniform(0, 1000) == priors.S
    @test Uniform(0, 0.1) == priors.B_gerda_all_pII
    α_min = maximum(-0.476981 ./ 0.0395483)
    @test Truncated(Normal(0, 1), α_min, Inf) == priors.αe_all
    prior_ω =
        Vector{Truncated{Normal{Float64},Continuous,Float64,Float64,Float64}}(undef, 1)
    prior_ω[1] = Truncated(Normal(1.328561380042463, 0.08067940552016985), 0, Inf)
    @test prior_ω == priors.ω.v
    prior_𝛥 =
        Vector{Truncated{Normal{Float64},Continuous,Float64,Float64,Float64}}(undef, 1)
    prior_𝛥[1] = Truncated(Normal(-0.33885135, 0.07638651), -Inf, Inf)
    @test prior_𝛥 == priors.𝛥.v
    # pretty names
    long_name = "GERDA part00 ANG4"
    @test "Energy Resolution " * L"(\omega)" * " " * long_name * " [keV]" ==
          pretty_names[:ω][1]
    @test "Energy Scale Bias " * L"(\Delta)" * " - " * long_name * " [keV]" ==
          pretty_names[:𝛥][1]
    @test L"\alpha_{\varepsilon} (" * "all)" == pretty_names[:αe_all]
    # nuisance_info
    @test [] == nuisance_info["α"]
    @test [] == nuisance_info["αr"]
    @test [] == nuisance_info["αb"]
    @test [] == nuisance_info["γ"]
    @test [] == nuisance_info["ε"]
    @test ["GERDA", "part00", "ANG4", 1.328561380042463, 0.08067940552016985, 0, Inf] ==
          nuisance_info["ω"][1]
    @test ["GERDA", "part00", "ANG4", -0.33885135, 0.07638651, -Inf, Inf] ==
          nuisance_info["𝛥"][1]
    @test ["combined", "", "", 0, 1, α_min, Inf] == nuisance_info["αe_all"][1]


    # correlated energy scale parameters 
    config["nuisance"]["energy_bias"]["correlated"] = true
    config["nuisance"]["energy_res"]["correlated"] = true
    partitions, fit_ranges = ZeroNuFit.Utils.get_partitions(config)
    events = ZeroNuFit.Utils.get_events(config["events"][1], partitions)
    part_event_index = ZeroNuFit.Utils.get_partition_event_index(events, partitions)
    settings = ZeroNuFit.Utils.get_settings(config)
    bkg_name =
        JSON.parsefile(config["partitions"][1])["fit_groups"]["all_phase_II"]["bkg_name"]
    priors = nothing
    pretty_names = nothing
    nuisance_info = nothing
    priors, pretty_names, nuisance_info =
        ZeroNuFit.Likelihood.build_prior(partitions, part_event_index, config, settings)

    # priors
    @test Uniform(0, 1000) == priors.S
    @test Uniform(0, 0.1) == priors.B_gerda_all_pII
    α_min = maximum(-0.476981 ./ 0.0395483)
    @test Truncated(Normal(0, 1), α_min, Inf) == priors.αe_all
    @test Truncated(Normal(0, 1), -Inf, Inf) == priors.αb_all
    αr_min = maximum(-1.328561380042463 ./ 0.08067940552016985)
    @test Truncated(Normal(0, 1), αr_min, Inf) == priors.αr_all
    # pretty names
    long_name = "GERDA part00 ANG4"
    @test L"\alpha_{\varepsilon} (" * "all)" == pretty_names[:αe_all]
    @test L"\alpha_{b} (" * "all)" == pretty_names[:αb_all]
    @test L"\alpha_{r} (" * "all)" == pretty_names[:αr_all]
    # nuisance_info
    @test [] == nuisance_info["α"]
    @test [] == nuisance_info["αr"]
    @test [] == nuisance_info["αb"]
    @test [] == nuisance_info["γ"]
    @test [] == nuisance_info["ε"]
    @test ["combined", "", "", 0, 1, -12.060720688373456, Inf] == nuisance_info["αe_all"][1]
    @test ["combined", "", "", 0, 1, -16.467168684210527, Inf] == nuisance_info["αr_all"][1]
    @test ["combined", "", "", 0, 1, -Inf, Inf] == nuisance_info["αb_all"][1]

    # no signal = no S prior
    config["bkg_only"] = true
    partitions, fit_ranges = ZeroNuFit.Utils.get_partitions(config)
    events = ZeroNuFit.Utils.get_events(config["events"][1], partitions)
    part_event_index = ZeroNuFit.Utils.get_partition_event_index(events, partitions)
    settings = ZeroNuFit.Utils.get_settings(config)
    bkg_name =
        JSON.parsefile(config["partitions"][1])["fit_groups"]["all_phase_II"]["bkg_name"]
    priors = nothing
    pretty_names = nothing
    nuisance_info = nothing
    priors, pretty_names, nuisance_info =
        ZeroNuFit.Likelihood.build_prior(partitions, part_event_index, config, settings)
    @test !haskey(priors, :S)
end
