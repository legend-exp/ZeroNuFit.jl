using Random, LinearAlgebra, Statistics, Distributions, StatsBase
using BAT, DensityInterface, IntervalSets

using TypedTables
using Plots, LaTeXStrings
using Cuba
using OrderedCollections

function get_bkg_pdf(
    bkg_shape::Symbol,
    evt_energy::Float64,
    p::NamedTuple,
    b_name::Symbol,
    fit_range,
)
    if (bkg_shape == :uniform)
        return norm_uniform(evt_energy, p, b_name, fit_range)
    elseif (bkg_shape == :linear)
        return norm_linear(evt_energy, p, b_name, fit_range)
    elseif (bkg_shape == :exponential)
        return norm_exponential(evt_energy, p, b_name, fit_range)
    else
        @error "bkg shape", bkg_shape, " is not yet implemented"
        exit(-1)
    end

end

function get_signal_pdf(evt_energy::Float64, Qbb::Float64, part_k::NamedTuple)
    signal_shape = part_k.signal_name
    bias = part_k.bias
    reso = part_k.width
    #@info part_k.experiment, signal_shape
    if (signal_shape == :gaussian)
        return pdf(Normal(Qbb - bias, reso), evt_energy)
    elseif (signal_shape == :gaussian_plus_lowEtail)
        return gaussian_plus_lowEtail(evt_energy, Qbb, bias, reso, part_k)
    else
        @error "signal shape ", signal_shape, " is not yet implemented"
        exit(-1)
    end

end


"""
    get_mu_b(deltaE, exposure, bkg)

Get the expected number of background counts in a partition
"""
function get_mu_b(deltaE, exposure, bkg_index)
    return deltaE * exposure * bkg_index
end

"""
    get_mu_s(deltaE, exposure, bkg)

Get the expected number of signal counts in a partition
"""
function get_mu_s(exposure, eff, signal)
    N_A = constants.N_A
    m_76 = constants.m_76
    sig_units = constants.sig_units
    return log(2) * N_A * exposure * (eff) * (signal * sig_units) / m_76
end


"""
    get_mu_s_b(p::NamedTuple,part_k::NamedTuple,idx_part_with_events::Int,settings::Dict,fit_range)

Get the expected number of signal and background counts in a partition
"""
function get_mu_s_b(
    p::NamedTuple,
    part_k::NamedTuple,
    idx_part_with_events::Int,
    settings::Dict,
    fit_range,
)

    deltaE = get_deltaE(fit_range)
    eff = get_efficiency(p, part_k, idx_part_with_events, settings)

    if (settings[:bkg_only] == false)
        model_s_k = get_mu_s(part_k.exposure, eff, p.S)
    else
        model_s_k = 0
    end

    b_name = part_k.bkg_name
    model_b_k = get_mu_b(deltaE, part_k.exposure, p[b_name])

    return model_s_k, model_b_k
end


"""
    build_likelihood_zero_obs_evts(part_k::NamedTuple, p::NamedTuple,settings::Dict,fit_range)

Function to calculate the partial likelihood for a partition with 0 events
"""
function build_likelihood_zero_obs_evts(
    part_k::NamedTuple,
    p::NamedTuple,
    settings::Dict,
    fit_range,
)

    ll_value = 0
    model_s_k, model_b_k = get_mu_s_b(p, part_k, 0, settings, fit_range)
    model_tot_k = model_b_k + model_s_k

    ll_value += -(model_tot_k + eps(model_tot_k))

    return ll_value
end


"""
    build_likelihood_per_partition(idx_k::Int, idx_part_with_events::Int,part_k::NamedTuple, events_k::Vector{Union{Float64}},p::NamedTuple,settings::Dict,bkg_shape::Symbol,fit_range)

Function which computes the partial likelihood for a single data partiton
free parameters: signal (S), background (B), energy bias (biask) and resolution per partition (resk)
"""
function build_likelihood_per_partition(
    idx_k::Int,
    idx_part_with_events::Int,
    part_k::NamedTuple,
    events_k::Vector{Union{Float64}},
    p::NamedTuple,
    settings::Dict,
    bkg_shape::Symbol,
    fit_range,
)
    Qbb = constants.Qbb

    ll_value = 0

    model_s_k, model_b_k = get_mu_s_b(p, part_k, idx_part_with_events, settings, fit_range)

    model_tot_k = model_b_k + model_s_k

    # constrain λ not to be negative
    if model_tot_k < 0
        λ = eps(0.0)
    else
        λ = model_tot_k + eps(model_tot_k)
    end

    ll_value += logpdf(Poisson(λ), length(events_k))

    for evt_energy in events_k

        term1 =
            model_b_k * get_bkg_pdf(bkg_shape, evt_energy, p, part_k.bkg_name, fit_range)

        if (settings[:bkg_only] == false)

            # get the correct reso and bias 
            reso, bias = get_energy_scale_pars(part_k, p, settings, idx_part_with_events)
            term2 = model_s_k * get_signal_pdf(evt_energy, Qbb, part_k)
        else
            term2 = 0
        end

        ll_value +=
            log((term1 + term2) + eps(term1 + term2)) - log(model_tot_k + eps(model_tot_k))

    end

    return ll_value
end



"""
    build_likelihood_looping_partitions(partitions::TypedTables.Table,events::Array{Vector{Float64}},part_event_index::Vector{Int},settings::Dict,sqrt_prior::Bool,s_max::Union{Float64,Nothing},fit_ranges;bkg_shape::Symbol=:uniform)

Function which creates the likelihood function for the fit (looping over partitions)

# Parameters
- `partitions::TypedTables.Table`: The partitions input table.
- `events::Array{Vector{Float64}}`: A list of events (=energies) in each partition.
- `part_event_index::Vector{Int}`: The index mapping events to the partitions.
- `settings::Dict`: A dictionary of settings containing configuration for the likelihood calculation.
- `sqrt_prior::Bool`: Whether to include the square root prior in the likelihood calculation. If `False`, a uniform prior is used.
- `s_max::Union{Float64, Nothing}`: A maximum value used for scaling the square root prior. If `Nothing`, no prior is applied.
- `fit_ranges`: The fitting ranges corresponding to the partitions.
- `bkg_shape::Symbol`: Specifies the background shape; default is `:uniform`.

Returns the likelihood function (a `DensityInterface.logfuncdensity` object).
"""
function build_likelihood_looping_partitions(
    partitions::TypedTables.Table,
    events::Array{Vector{Float64}},
    part_event_index::Vector{Int},
    settings::Dict,
    sqrt_prior::Bool,
    s_max::Union{Float64,Nothing},
    fit_ranges;
    bkg_shape::Symbol = :uniform,
)
    @debug part_event_index
    return DensityInterface.logfuncdensity(
        function (p::NamedTuple)
            total_ll = 0.0

            for (idx_k, part_k) in enumerate(partitions)

                if part_event_index[idx_k] != 0
                    idx_k_with_events = part_event_index[idx_k]
                    total_ll += build_likelihood_per_partition(
                        idx_k,
                        part_event_index[idx_k],
                        part_k,
                        events[idx_k],
                        p,
                        settings,
                        bkg_shape,
                        fit_ranges[part_k.fit_group],
                    )
                else
                    # no events are there for a given partition
                    total_ll += build_likelihood_zero_obs_evts(
                        part_k,
                        p,
                        settings,
                        fit_ranges[part_k.fit_group],
                    )
                end
            end

            ## trick to include thes sqrt prior 
            if (sqrt_prior)
                total_ll += -log(2) - 0.5 * log(s_max) - 0.5 * log(p.S + eps(p.S))
            end

            return total_ll
        end,
    )
end




"""
    generate_data(samples::BAT.DensitySampleVector,partitions::TypedTables.Table,part_event_index::Vector{Int},settings::Dict,fit_ranges;best_fit::Bool=false,seed=nothing,bkg_only=false)

Generates data from a posterior distribution.
This is based on the posterior predictive distributions. 
Given a model with some parameters `theta_i`, the posterior predictive distribution,
or the distribution of data generated according to the posterior distribution of theta
and the likelihood is:

```math
p(y|D) =int p(y|theta)p(theta|D)dtheta
```

Or in terms of sampling we first draw samples of `theta` from the posterior and then generate,
datasets based on the likelihood.
We also give the options to fix the posterior distribution to the best fit,
which is equivalent to the standard sampling methods.

Parameters
----------
    - samples::DensitySamplesVector the samples of a past fit or a NamedTuple of best fit
    - partitions::Table of the partition info
    - part_event_index: index for the parameters for partitions with events

Keyword arguments
-----------------
    - best_fit::Bool where to fix the paramaters to the best fit
    - nuis_prior::Bool whether only statistical parameters were included in the posterior
    - bkg_only::Bool where the fit was without signal,
    - seed::Int random seed

Returns
-------
    OrderedDict of the data
"""
function generate_data(
    samples::BAT.DensitySampleVector,
    partitions::TypedTables.Table,
    part_event_index::Vector{Int},
    settings::Dict,
    fit_ranges;
    best_fit::Bool = false,
    seed = nothing,
    bkg_only = false,
)
    Qbb = constants.Qbb

    # seed the seed
    output = OrderedDict("events" => [])
    if (seed == nothing)
        Random.seed!(Int(round(10000 * (rand()))))
    else
        Random.seed!(seed)
    end


    if (samples isa NamedTuple)
        p = samples
    else
        distribution = Categorical(samples.weight / sum(samples.weight)) # Generate a random index based on the weights
        random_index = rand(distribution)
        p = samples.v[random_index]
    end

    # create the array to fill
    events = []

    for (idx_k, part_k) in enumerate(partitions)

        b_name = part_k.bkg_name
        idx_part_with_events = part_event_index[idx_k]
        model_s_k, model_b_k = get_mu_s_b(
            p,
            part_k,
            idx_part_with_events,
            settings,
            fit_ranges[part_k.fit_group],
        )

        n_s = rand(Poisson(model_s_k))
        n_b = rand(Poisson(model_b_k))
        events = generate_disjoint_uniform_samples(n_b, fit_ranges[part_k.fit_group])
        if (bkg_only == false)
            for i = 1:n_s

                reso, bias =
                    get_energy_scale_pars(part_k, p, settings, idx_part_with_events)

                append!(events, rand(Normal(Qbb - bias, reso)))


            end
        end
        times = rand(Uniform(part_k.start_ts, part_k.end_ts), length(events))
        times = [Int(round(t)) for t in times]

        if (length(times) > 0)
            for (t, e) in zip(times, events)
                append!(
                    output["events"],
                    [
                        OrderedDict(
                            "timestamp" => t,
                            "experiment" => part_k.experiment,
                            "detector" => part_k.detector,
                            "energy" => e,
                        ),
                    ],
                )
            end

        end

    end
    display(output["events"])

    return output

end



"""
    get_signal_bkg_priors(config)

Defines specific priors for signal and background contributions

Parameters
----------
    - config: the Dict of the fit config
"""
function get_signal_bkg_priors(config)

    uppS = config["signal"]["upper_bound"]
    uppB = config["bkg"]["upper_bound"]

    if config["signal"]["prior"] == "uniform" || config["signal"]["prior"] == "sqrt"
        distrS = 0 .. uppS
    elseif config["signal"]["prior"] == "loguniform"
        distrS = LogUniform(0.01, uppS)
    else
        @error "Distribution", config["signal"]["prior"], " is not yet defined"
        exit(-1)
    end

    if config["bkg"]["prior"] == "uniform"
        distrB = 0 .. uppB
    else
        @error "Distribution", config["bkg"]["prior"], " is not yet defined"
        exit(-1)
    end

    return distrS, distrB
end


"""
    build_prior(partitions,part_event_index,config,settings::Dict;hierachical=false,hierachical_mode=nothing,hierachical_range=nothing,bkg_shape=:uniform,shape_pars=nothing)

Builds the priors for use in the fit

Parameters
----------
    - partitions:Table of the partition info
    - config: the Dict of the fit config
    - nuis_prior; true if we want to include priors for nuisance parameters (bias, res, eff)
"""
function build_prior(
    partitions,
    part_event_index,
    config,
    settings::Dict;
    hierachical = false,
    hierachical_mode = nothing,
    hierachical_range = nothing,
    bkg_shape = :uniform,
    shape_pars = nothing,
)

    # bkg indexs
    list_names = partitions.bkg_name
    unique_list = unique(list_names)
    bkg_names = [Symbol(name) for name in unique_list]

    distrS, distrB = get_signal_bkg_priors(config)
    distrB_multi = OrderedDict(Symbol(bkg_name) => distrB for bkg_name in bkg_names)

    pretty_names = OrderedDict(
        :S => string("S [") * L"10^{-27}" * string("yr") * L"^{-1}" * string("]"),
        :α => L"\alpha_{\varepsilon}",
        :αr => L"\alpha_{r}",
        :αb => L"\alpha_{b}",
        :γ => [],
        :ε => [],
        :ω => [],
        :𝛥 => [],
    )

    for key in keys(distrB_multi)
        pretty_names[key] = string(key) * " [cts/keV/kg/yr]"
    end


    # create priors one by one
    ### SIGNAL PRIOR

    priors = OrderedDict()
    if (settings[:bkg_only] == false)
        priors[:S] = distrS
        @info "entered to add S prior"
    end

    # dictionary with info on the prior parameters
    nuisance_info = OrderedDict(
        "α" => [],
        "αr" => [],
        "αb" => [],
        "γ" => [],
        "ε" => [],
        "ω" => [],
        "𝛥" => [],
    )

    ### EFFICIENCY prior

    if (settings[:eff_fixed] == false && settings[:eff_correlated] == true)

        all_eff_tot = partitions.eff_tot
        all_eff_tot_sigma = partitions.eff_tot_sigma
        ratio = -all_eff_tot ./ all_eff_tot_sigma
        α_min = maximum(ratio)

        list_names = partitions.eff_name
        unique_list = unique(list_names)
        for name in unique_list
            priors[Symbol(name)] = Truncated(Normal(0, 1), α_min, Inf)
            pretty_names[Symbol(name)] =
                L"\alpha_{\varepsilon} (" * split(String(name), "_")[2] * ")"
            nuisance_info[string(name)] = [["combined", "", "", 0, 1, α_min, Inf]]
        end
    end
    if (
        settings[:eff_fixed] == false &&
        settings[:eff_correlated] == false
    )
        eff = Vector{Truncated{Normal{Float64},Continuous,Float64,Float64,Float64}}(
            undef,
            maximum(part_event_index),
        )

        for (idx, part) in enumerate(partitions)

            if (part_event_index[idx] != 0)
                i_new = part_event_index[idx]

                eff[i_new] = Truncated(Normal(part.eff_tot, part.eff_tot_sigma), 0, 1)
                long_name =
                    string(part.experiment) *
                    " " *
                    string(part.part_name) *
                    " " *
                    part.detector
                append!(
                    pretty_names[:ε],
                    ["Efficiency " * L"(\varepsilon)" * " - " * long_name * ""],
                )
                append!(
                    nuisance_info["ε"],
                    [[
                        string(part.experiment),
                        string(part.part_name),
                        part.detector,
                        part.eff_tot,
                        part.eff_tot_sigma,
                        0,
                        1,
                    ]],
                )
            end
        end
        priors[:ε] = eff
    end

    ### ENERGY BIAS prior

    if (
        settings[:energy_bias_fixed] == false &&
        settings[:energy_bias_correlated] == true
    )

        list_names = partitions.energy_bias_name
        unique_list = unique(list_names)
        for name in unique_list
            priors[Symbol(name)] = Truncated(Normal(0, 1), -Inf, Inf)
            pretty_names[Symbol(name)] = L"\alpha_{b} (" * split(String(name), "_")[2] * ")"
            nuisance_info[string(name)] =
                [["combined", "", "", part.detector, 0, 1, -Inf, Inf]]

        end
    end
    if (
        settings[:energy_bias_fixed] == false &&
        settings[:energy_bias_correlated] == false
    )
        bias = Vector{Truncated{Normal{Float64},Continuous,Float64,Float64,Float64}}(
            undef,
            maximum(part_event_index),
        )

        for (idx, part) in enumerate(partitions)

            if (part_event_index[idx] != 0)
                i_new = part_event_index[idx]

                long_name =
                    string(part.experiment) *
                    " " *
                    string(part.part_name) *
                    " " *
                    part.detector
                if part.signal_name == :gaussian
                    bias[i_new] = Truncated(Normal(part.bias, part.bias_sigma), -Inf, Inf)
                    append!(
                        nuisance_info["𝛥"],
                        [[
                            string(part.experiment),
                            string(part.part_name),
                            part.detector,
                            part.bias,
                            part.bias_sigma,
                            -Inf,
                            Inf,
                        ]],
                    )

                elseif part.signal_name == :gaussian_plus_lowEtail
                    # let's define some intervals in +-5σ (always with res>0)
                    bias_min = part.bias - 5 * part.bias_sigma
                    bias_max = part.bias + 5 * part.bias_sigma
                    bias[i_new] =
                        Truncated(Normal(part.bias, part.bias_sigma), bias_min, bias_max)
                    append!(
                        nuisance_info["𝛥"],
                        [[
                            string(part.experiment),
                            string(part.part_name),
                            part.detector,
                            part.bias,
                            part.bias_sigma,
                            bias_min,
                            bias_max,
                        ]],
                    )

                else
                    @info "There is no specific name for the $(part.signal_name) peak shape, exit here"
                    exit()
                end
                append!(
                    pretty_names[:𝛥],
                    ["Energy Scale Bias " * L"(\Delta)" * " - " * long_name * " [keV]"],
                )
            end
        end
        priors[:𝛥] = bias
    end

    ### ENERGY RESOLUTION prior
    if (
        settings[:energy_res_fixed] == false &&
        settings[:energy_res_correlated] == true
    )
        all_width = partitions.width
        all_width_sigma = partitions.width_sigma
        ratio = -all_width ./ all_width_sigma
        αr_min = maximum(ratio)

        list_names = partitions.energy_reso_name
        unique_list = unique(list_names)
        for name in unique_list
            priors[Symbol(name)] = Truncated(Normal(0, 1), αr_min, Inf)
            pretty_names[Symbol(name)] = L"\alpha_{r} (" * split(String(name), "_")[2] * ")"
            nuisance_info[string(name)] = [["combined", "", "", 0, 1, αr_min, Inf]]

        end

    
    end
    if (
        settings[:energy_res_fixed] == false &&
        settings[:energy_res_correlated] == false
    )
        res = Vector{Truncated{Normal{Float64},Continuous,Float64,Float64,Float64}}(
            undef,
            maximum(part_event_index),
        )

        for (idx, part) in enumerate(partitions)

            if (part_event_index[idx] != 0)
                i_new = part_event_index[idx]

                long_name =
                    string(part.experiment) *
                    " " *
                    string(part.part_name) *
                    " " *
                    part.detector
                if part.signal_name == :gaussian
                    res[i_new] = Truncated(Normal(part.width, part.width_sigma), 0, Inf)
                    append!(
                        pretty_names[:ω],
                        ["Energy Resolution " * L"(\omega)" * " " * long_name * " [keV]"],
                    )
                    append!(
                        nuisance_info["ω"],
                        [[
                            string(part.experiment),
                            string(part.part_name),
                            part.detector,
                            part.width,
                            part.width_sigma,
                            0,
                            Inf,
                        ]],
                    )

                elseif part.signal_name == :gaussian_plus_lowEtail
                    # let's define some intervals in +-5σ (always with res>0)
                    res_min = part.width - 5 * part.width_sigma
                    res_min = res_min < 0 ? 0 : res_min
                    res_max = part.width + 5 * part.width_sigma
                    res[i_new] =
                        Truncated(Normal(part.width, part.width_sigma), res_min, res_max)
                    append!(
                        pretty_names[:ω],
                        [
                            "Resolution fractional uncertainty " *
                            L"(\omega)" *
                            " " *
                            long_name *
                            " [keV]",
                        ],
                    )
                    append!(
                        nuisance_info["ω"],
                        [[
                            string(part.experiment),
                            string(part.part_name),
                            part.detector,
                            part.width,
                            part.width_sigma,
                            res_min,
                            res_max,
                        ]],
                    )

                else
                    @info "There is no specific name for the $(part.signal_name) peak shape, exit here"
                    exit()
                end
            end
        end
        priors[:ω] = res

    end

    ## bkg shape priors
    if shape_pars != nothing

        for par in keys(shape_pars)
            name = par
            prior = shape_pars[par]

            for bkg_name in bkg_names
                priors[Symbol(string(bkg_name) * "_" * name)] = Uniform(prior[1], prior[2])
                pretty_names[Symbol(string(bkg_name) * "_" * name)] =
                    string(bkg_name) * "_" * name
            end

        end

    end


    ## BKG prior
    if (hierachical == false)
        for (key, item) in distrB_multi
            priors[key] = item
        end

        return distprod(; priors...), pretty_names, nuisance_info
    else

        if (hierachical_mode == "lognormal")
            hd = BAT.HierarchicalDistribution(
                v -> begin
                    dict = (;
                        (
                            key => LogNormal(log(v.B) - 0.5 * v.σB * v.σB, v.σB) for
                            key in keys(distrB_multi)
                        )...
                    )
                    return distprod(; dict..., priors...)
                end,
                distprod(
                    B = 0 .. 1,
                    σB = Uniform(hierachical_range[1], hierachical_range[2]),
                ),
            )
        elseif (hierachical_mode == "normal")
            hd = BAT.HierarchicalDistribution(
                v -> begin
                    dict = (;
                        (
                            key => Normal(log(v.B) - 0.5 * v.σB * v.σB, v.σB) for
                            key in keys(distrB_multi)
                        )...
                    )
                    return distprod(; dict..., priors...)
                end,
                distprod(
                    B = 0 .. 1,
                    σB = Uniform(hierachical_range[1], hierachical_range[2]),
                ),
            )
        else
            @error "hierachical (correlated) bkg mode $hierachical_mode is not know"
            exit(-1)
        end

        pretty_names[:B] = "B [cts/keV/kg/yr]"
        pretty_names[:σB] = L"\sigma_B" * string("[cts/keV/kg/yr]")

        return hd, pretty_names, nuisance_info
    end
end
