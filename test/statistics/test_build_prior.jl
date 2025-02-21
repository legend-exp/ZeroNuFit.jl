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

Base.exit(code::Int) = throw(ArgumentError("exit code $code"))

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

    # eff
    config["nuisance"]["efficiency"]["fixed"] = false
    config["nuisance"]["efficiency"]["correlated"] = false
    partitions, fit_ranges = ZeroNuFit.Utils.get_partitions(config)
    events = ZeroNuFit.Utils.get_events(config["events"][1], partitions)
    part_event_index = ZeroNuFit.Utils.get_partition_event_index(events, partitions)
    settings = ZeroNuFit.Utils.get_settings(config)
    priors = nothing
    pretty_names = nothing
    nuisance_info = nothing
    priors, pretty_names, nuisance_info =
        ZeroNuFit.Likelihood.build_prior(partitions, part_event_index, config, settings)
    # priors
    prior_ε =
        Vector{Truncated{Normal{Float64},Continuous,Float64,Float64,Float64}}(undef, 1)
    prior_ε[1] = Truncated(Normal(0.476981, 0.0395483), 0, 1)
    @test prior_ε == priors.ε.v
    # pretty names
    @test "Efficiency " * L"(\varepsilon)" * " - " * long_name * "" == pretty_names[:ε][1]
    # nuisance info
    @test ["GERDA", "part00", "ANG4", 0.476981, 0.0395483, 0, 1] == nuisance_info["ε"][1]


    # bkg shapes parameters (eg linear/exponential)
    config["bkg"] = Dict(
        "correlated" => Dict("range" => "none", "mode" => "none"),
        "prior" => "uniform",
        "upper_bound" => 0.1,
        "shape" => Dict("name" => "linear", "pars" => Dict("slope" => [-10, 10])),
    )
    partitions, fit_ranges = ZeroNuFit.Utils.get_partitions(config)
    events = ZeroNuFit.Utils.get_events(config["events"][1], partitions)
    part_event_index = ZeroNuFit.Utils.get_partition_event_index(events, partitions)
    settings = ZeroNuFit.Utils.get_settings(config)
    bkg_shape, bkg_shape_pars = ZeroNuFit.Utils.get_bkg_info(config)
    priors = nothing
    pretty_names = nothing
    nuisance_info = nothing
    priors, pretty_names, nuisance_info = ZeroNuFit.Likelihood.build_prior(
        partitions,
        part_event_index,
        config,
        settings,
        shape_pars = bkg_shape_pars,
    )
    # priors
    @test Uniform(-10, 10) == priors.B_gerda_all_pII_slope
    # pretty names
    long_name = "GERDA part00 ANG4"
    @test bkg_name * "_slope" == pretty_names[:B_gerda_all_pII_slope]

    # gaussian + low E tail (MJD-like event)
    config["bkg"] = Dict(
        "correlated" => Dict("range" => "none", "mode" => "none"),
        "prior" => "uniform",
        "upper_bound" => 0.1,
    )
    config["nuisance"]["energy_bias"]["fixed"] = false
    config["nuisance"]["energy_bias"]["correlated"] = false
    config["nuisance"]["energy_res"]["fixed"] = false
    config["nuisance"]["energy_res"]["correlated"] = false
    config["events"] = [joinpath(present_dir, "../inputs/events_test_2.json")]
    config["partitions"] = [joinpath(present_dir, "../inputs/partitions_test_2.json")]
    partitions, fit_ranges = ZeroNuFit.Utils.get_partitions(config)
    events = ZeroNuFit.Utils.get_events(config["events"][1], partitions)
    part_event_index = ZeroNuFit.Utils.get_partition_event_index(events, partitions)
    settings = ZeroNuFit.Utils.get_settings(config)
    priors = nothing
    pretty_names = nothing
    nuisance_info = nothing
    priors, pretty_names, nuisance_info =
        ZeroNuFit.Likelihood.build_prior(partitions, part_event_index, config, settings)
    # priors
    prior_𝛥 =
        Vector{Truncated{Normal{Float64},Continuous,Float64,Float64,Float64}}(undef, 1)
    prior_𝛥[1] = Truncated(Normal(0.02, 0.08), 0.02 - 5 * 0.08, 0.02 + 5 * 0.08)
    @test prior_𝛥 == priors.𝛥.v
    prior_ω =
        Vector{Truncated{Normal{Float64},Continuous,Float64,Float64,Float64}}(undef, 1)
    prior_ω[1] = Truncated(Normal(1.05, 0.040), 1.05 - 5 * 0.040, 1.05 + 5 * 0.040)
    @test prior_ω == priors.ω.v
    # nuisance info
    @test ["MJD", "ds0", "DS0", 0.02, 0.08, 0.02 - 5 * 0.08, 0.02 + 5 * 0.08] ==
          nuisance_info["𝛥"][1]
    @test ["MJD", "ds0", "DS0", 1.05, 0.04, 1.05 - 5 * 0.040, 1.05 + 5 * 0.040] ==
          nuisance_info["ω"][1]

    # same, but with a not-existing signal shape
    config["partitions"] = [joinpath(present_dir, "../inputs/partitions_test_4.json")]
    partitions, fit_ranges = ZeroNuFit.Utils.get_partitions(config)
    events = ZeroNuFit.Utils.get_events(config["events"][1], partitions)
    part_event_index = ZeroNuFit.Utils.get_partition_event_index(events, partitions)
    settings = ZeroNuFit.Utils.get_settings(config)
    @test_throws ArgumentError ZeroNuFit.Likelihood.build_prior(
        partitions,
        part_event_index,
        config,
        settings,
    )

    # hierarchical (back to GERDA-like event)
    config["bkg"]["correlated"] = Dict("mode" => "lognormal", "range" => [0, 1])
    config["events"] = [joinpath(present_dir, "../inputs/events_test.json")]
    config["partitions"] = [joinpath(present_dir, "../inputs/partitions_test.json")]
    partitions, fit_ranges = ZeroNuFit.Utils.get_partitions(config)
    events = ZeroNuFit.Utils.get_events(config["events"][1], partitions)
    part_event_index = ZeroNuFit.Utils.get_partition_event_index(events, partitions)
    settings = ZeroNuFit.Utils.get_settings(config)
    corr, hier_mode, hier_range = ZeroNuFit.Utils.get_corr_info(config)
    priors = nothing
    pretty_names = nothing
    nuisance_info = nothing
    priors, pretty_names, nuisance_info = ZeroNuFit.Likelihood.build_prior(
        partitions,
        part_event_index,
        config,
        settings,
        hierachical = corr,
        hierachical_mode = hier_mode,
        hierachical_range = hier_range,
    )
    # pretty names
    @test "B [cts/keV/kg/yr]" == pretty_names[:B]
    @test L"\sigma_B" * string("[cts/keV/kg/yr]") == pretty_names[:σB]

    # unknown hierarchical mode
    config["bkg"]["correlated"] = Dict("mode" => "fancynormal", "range" => [0, 1])
    partitions, fit_ranges = ZeroNuFit.Utils.get_partitions(config)
    events = ZeroNuFit.Utils.get_events(config["events"][1], partitions)
    part_event_index = ZeroNuFit.Utils.get_partition_event_index(events, partitions)
    settings = ZeroNuFit.Utils.get_settings(config)
    corr, hier_mode, hier_range = ZeroNuFit.Utils.get_corr_info(config)
    @test_throws ArgumentError ZeroNuFit.Likelihood.build_prior(
        partitions,
        part_event_index,
        config,
        settings,
        hierachical = corr,
        hierachical_mode = hier_mode,
        hierachical_range = hier_range,
    )
end
