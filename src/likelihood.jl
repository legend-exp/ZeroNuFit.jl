### likelihood.jl
#
# Authors: Sofia Calgaro, Toby Dixon
# 
###
module Likelihood
using Random, LinearAlgebra, Statistics, Distributions, StatsBase
using BAT, DensityInterface, IntervalSets

using TypedTables
using Plots, LaTeXStrings
using Cuba
using OrderedCollections
using SpecialFunctions

using ZeroNuFit

"""
    get_stat_blocks(partitions,events::Array{Vector{Float64}},part_event_index,fit_ranges;config,bkg_only::Bool)

Function to retrieve useful pieces (prior, likelihood, posterior), also in saving values.

### Arguments
- `partitions`: table of partitions.
- `events::Array{Vector{Float64}}`: list of events (=energies) in each partition.
- `part_event_index::Vector{Int}`: index mapping events to the partitions.
- `fit_ranges`: dictionary of energy ranges considered for the analysis.
- `config`: input dictionary.
- `bkg_only::Bool`: True if we are using a model with background only.
"""
function get_stat_blocks(
    partitions,
    events::Array{Vector{Float64}},
    part_event_index,
    fit_ranges;
    config,
    bkg_only::Bool,
)
    settings = ZeroNuFit.Utils.get_settings(config)

    corr, hier_mode, hier_range = ZeroNuFit.Utils.get_corr_info(config)

    bkg_shape, bkg_shape_pars = ZeroNuFit.Utils.get_bkg_info(config)

    prior, par_names, nuisance_info = build_prior(
        partitions,
        part_event_index,
        config,
        settings,
        hierachical = corr,
        hierachical_mode = hier_mode,
        hierachical_range = hier_range,
        bkg_shape = bkg_shape,
        shape_pars = bkg_shape_pars,
    )
    @info "using a ", bkg_shape, " bkg with ", bkg_shape_pars, " parameters"
    @info "built prior"

    sqrt_prior, s_max = get_signal_prior_info(bkg_only, config)

    likelihood = build_likelihood_looping_partitions(
        partitions,
        events,
        part_event_index,
        settings,
        sqrt_prior,
        s_max,
        fit_ranges,
        config["bkg"]["units"],
        bkg_shape = bkg_shape,
    )
    @info "built likelihood"

    posterior = PosteriorMeasure(likelihood, prior)
    @info "got posterior"

    return prior, likelihood, posterior, par_names, nuisance_info
end


"""
    run_fit_over_partitions(partitions,events::Array{Vector{Float64}},part_event_index::Vector{Int}, config,fit_ranges)

Function to run the fit looping over partitions.

### Arguments
- `partitions`: table of partitions.
- `events::Array{Vector{Float64}}`: list of events (=energies) in each partition.
- `part_event_index::Vector{Int}`: index mapping events to the partitions.
- `config`: input dictionary.
- `fit_ranges`: dictionary of energy ranges considered for the analysis.
"""
function run_fit_over_partitions(
    partitions,
    events::Array{Vector{Float64}},
    part_event_index::Vector{Int},
    config,
    fit_ranges,
)
    bkg_only = config["bkg_only"]
    prior, likelihood, posterior, par_names, nuisance_info = get_stat_blocks(
        partitions,
        events,
        part_event_index,
        fit_ranges,
        config = config,
        bkg_only = bkg_only,
    )

    Ns = Int(config["bat_fit"]["nsteps"])
    Nc = Int(config["bat_fit"]["nchains"])
    return bat_sample(
        posterior,
        MCMCSampling(mcalg = MetropolisHastings(), nsteps = Ns, nchains = Nc),
    ).result,
    prior,
    par_names
end

"""
    get_signal_prior_info(bkg_only::Bool,config)

Function that retrieves signal prior information.

### Arguments
- `bkg_only::Bool`: True if we are using a model with background only.
- `config`: input dictionary.
"""
function get_signal_prior_info(bkg_only::Bool, config)
    sqrt_prior = false
    s_max = nothing
    if bkg_only == false
        if (config["signal"]["prior"] == "sqrt")
            sqrt_prior = true
            s_max = Float64(config["signal"]["upper_bound"])
        end
    end
    return sqrt_prior, s_max
end


"""
    norm_uniform(x::Real,p::NamedTuple,fit_range)

Normalised flat function defined by 1/norm.

### Arguments
- `x::Real`: the x value to evaluate at.
- `p::NamedTuple`: collection of key-value pairs where each key corresponds to a model parameter.
- `fit_range`: array of arrays, defining the allowed energy ranges; e.g. `fit_range= [[1930,1950], [1970,1990], [2000,2050]]`.
"""
function norm_uniform(x::Real, p::NamedTuple, fit_range)
    range_l, range_h = ZeroNuFit.Utils.get_range(fit_range)
    center = range_l[1]

    norm = sum(range_h .- range_l)
    return 1 / norm
end


"""
    norm_linear(x::Float64,p::NamedTuple,b_name::Symbol,fit_range)

Normalised linear function defined by (1+slope*(x-center)/net_width)/norm.

### Arguments
- `x::Real`: the x value to evaluate at.
- `p::NamedTuple`: collection of key-value pairs where each key corresponds to a model parameter.
- `b_name::Symbol`: name of the background index.
- `fit_range`: array of arrays, defining the allowed energy ranges; e.g. `fit_range= [[1930,1950], [1970,1990], [2000,2050]]`.
"""
function norm_linear(x::Float64, p::NamedTuple, b_name::Symbol, fit_range)
    range_l, range_h = ZeroNuFit.Utils.get_range(fit_range)
    center = range_l[1]

    sum_range = sum(range_h .- range_l)
    sum_range_sq = sum(range_h .^ 2 .- range_l .^ 2)
    slope = p[Symbol(string(b_name) * "_slope")]

    delta = range_h[end] - range_l[1]
    norm = sum_range * (1 - slope * center / delta) + slope * sum_range_sq / (2 * delta)

    return (1 + slope * (x - center) / delta) / norm
end


"""
    exp_stable(x::Float64)

Exponential function, using Taylor expansion series for `abs(x) < 1E-6`.

### Arguments
- `x::Real`: the x value to evaluate at.
"""
function exp_stable(x::Float64)
    if (abs(x) < 1E-6)
        return 1 + x + x^2 / 2 + x^3 / 6
    else
        return exp(x)
    end
end


"""
    norm_exponential(x::Float64,p::NamedTuple,b_name::Symbol,fit_range)

Normalised exponential function defined by exp_stable((x-center)*Rt)/norm.

### Arguments
- `x::Real`: the x value to evaluate at.
- `p::NamedTuple`: collection of key-value pairs where each key corresponds to a model parameter.
- `b_name::Symbol`: name of the background index.
- `fit_range`: array of arrays, defining the allowed energy ranges; e.g. `fit_range= [[1930,1950], [1970,1990], [2000,2050]]`.
"""
function norm_exponential(x::Float64, p::NamedTuple, b_name::Symbol, fit_range)
    range_l, range_h = ZeroNuFit.Utils.get_range(fit_range)
    center = range_l[1]

    centers = fill(center, length(range_l))
    R = p[Symbol(string(b_name) * "_slope")]
    delta = range_h[end] - range_l[1]
    Rt = R / delta

    if (abs(Rt) > 1E-6)
        norm =
            (
                -sum(exp_stable.((range_l - centers) * Rt)) +
                sum(exp_stable.((range_h - centers) * Rt))
            ) / Rt
    else
        norm = sum(range_h .- range_l)
    end

    return exp_stable((x - center) * Rt) / norm

end


"""
    gaussian_plus_lowEtail(evt_energy::Float64,Qbb::Float64,bias::Float64,reso::Float64,part_k::NamedTuple)

Signal model based on the peak shape used for the MJD analysis. The peak shape was derived from considerations made in [S. I. Alvis et al., Phys. Rev. C 100, 025501 (2019)].

### Arguments
- `evt_energy::Float64`: energy of the event at which we want to compute the function.
- `Qbb::Float64`: centroid of the Gaussian.
- `bias::Float64`: energy bias associated to the energy event.
- `reso::Float64`: energy resolution associated to the energy event.
- `part_k::NamedTuple`: Table of specifications for a given partition k.
"""
function gaussian_plus_lowEtail(
    evt_energy::Float64,
    Qbb::Float64,
    bias::Float64,
    reso::Float64,
    part_k::NamedTuple,
)
    γ = reso
    # following params are ALWAYS fixed
    f = part_k.frac
    τ = part_k.tau
    σ = part_k.sigma

    term1 = (1 - f) * pdf(Normal(Qbb - bias, γ * σ), evt_energy)

    term2 =
        f / (2 * γ * τ) *
        exp(((γ * σ)^2) / (2 * (γ * τ)^2) + (evt_energy - (Qbb - bias)) / (γ * τ))
    term2 =
        term2 * erfc(σ / (sqrt(2) * τ) + (evt_energy - (Qbb - bias)) / (sqrt(2) * γ * σ))

    return term1 + term2
end


"""
    get_bkg_pdf(bkg_shape::Symbol,evt_energy::Float64,p::NamedTuple,b_name::Symbol,fit_range)

Returns the background modeling function.

### Arguments
- `bkg_shape::Symbol`: Specifies the background shape; default is `:uniform`.
- `evt_energy::Float64`: energy of the event at which we want to compute the function.
- `p::NamedTuple`: collection of key-value pairs where each key corresponds to a model parameter.
- `b_name::Symbol`: name of the background index.
- `fit_range`: array of arrays, defining the allowed energy ranges; e.g. `fit_range= [[1930,1950], [1970,1990], [2000,2050]]`.
"""
function get_bkg_pdf(
    bkg_shape::Symbol,
    evt_energy::Float64,
    p::NamedTuple,
    b_name::Symbol,
    fit_range,
)
    if (bkg_shape == :uniform)
        return norm_uniform(evt_energy, p, fit_range)
    elseif (bkg_shape == :linear)
        return norm_linear(evt_energy, p, b_name, fit_range)
    elseif (bkg_shape == :exponential)
        return norm_exponential(evt_energy, p, b_name, fit_range)
    else
        @error "bkg shape", bkg_shape, " is not yet implemented"
        exit(-1)
    end

end


"""
    get_signal_pdf(evt_energy::Float64, Qbb::Float64, part_k::NamedTuple)

Returns the signal modeling function.

### Arguments
- `evt_energy::Float64`: energy of the event at which we want to compute the function.
- `Qbb::Float64`: centroid of the Gaussian.
- `part_k::NamedTuple`: Table of specifications for a given partition k.
"""
function get_signal_pdf(evt_energy::Float64, Qbb::Float64, part_k::NamedTuple)
    signal_shape = part_k.signal_name
    bias = part_k.bias
    reso = part_k.width

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
    get_mu_b(deltaE, exposure, bkg_index, reso, bkg_units::String)

Get the expected number of background counts.

### Arguments
- `deltaE`: net width of the fit range.
- `exposure`: exposure (mass x time).
- `bkg_index`: background index value.
- `bkg_units::String`: Specifies the units for the background index; available options are `"ckky"` (=counts/keV/kg/yr) or `"cFty"` (=counts/FWHM/t/yr).
"""
function get_mu_b(deltaE, exposure, bkg_index, reso, bkg_units::String)
    if bkg_units == "ckky"
        return deltaE * exposure * bkg_index
    end
    if bkg_units == "cFty"
        fwhm = reso * 2.355
        return deltaE * exposure * (bkg_index / fwhm / 1000)
    end
end

"""
    get_mu_s(exposure, eff, signal)

Get the expected number of signal counts.

### Arguments
- `exposure`: exposure (mass x time).
- `eff`: total signal efficiency.
- `signal`: number of signal counts
"""
function get_mu_s(exposure, eff, signal)
    N_A = ZeroNuFit.Constants.N_A
    m_76 = ZeroNuFit.Constants.m_76
    sig_units = ZeroNuFit.Constants.sig_units
    return log(2) * N_A * exposure * (eff) * (signal * sig_units) / m_76
end


"""
    get_mu_s_b(p::NamedTuple,part_k::NamedTuple,idx_part_with_events::Int,settings::Dict,fit_range,bkg_units::String)

Get the expected number of signal and background counts.

### Arguments
- `p::NamedTuple`: collection of key-value pairs where each key corresponds to a model parameter.
- `part_k::NamedTuple`: Table of specifications for a given partition k.
- `idx_part_with_events::Int`: index of the partition with the event.
- `settings::Dict`: dictionary of settings containing configuration for the likelihood calculation.
- `fit_range`: array of arrays, defining the allowed energy ranges; e.g. `fit_range= [[1930,1950], [1970,1990], [2000,2050]]`.
- `bkg_units::String`: Specifies the units for the background index; available options are `"ckky"` (=counts/keV/kg/yr) or `"cFty"` (=counts/FWHM/t/yr).
"""
function get_mu_s_b(
    p::NamedTuple,
    part_k::NamedTuple,
    idx_part_with_events::Int,
    settings::Dict,
    fit_range,
    bkg_units::String,
)

    deltaE = ZeroNuFit.Utils.get_deltaE(fit_range)
    eff = ZeroNuFit.Utils.get_efficiency(p, part_k, idx_part_with_events, settings)

    if (settings[:bkg_only] == false)
        model_s_k = get_mu_s(part_k.exposure, eff, p.S)
    else
        model_s_k = 0
    end

    b_name = part_k.bkg_name
    reso, _ =
        ZeroNuFit.Utils.get_energy_scale_pars(part_k, p, settings, idx_part_with_events)
    model_b_k = get_mu_b(deltaE, part_k.exposure, p[b_name], reso, bkg_units)

    return model_s_k, model_b_k
end


"""
    build_likelihood_zero_obs_evts(part_k::NamedTuple, p::NamedTuple,settings::Dict,fit_range,bkg_units::String)

Function to calculate the likelihood for a single data partition k with 0 events.

### Arguments
- `part_k::NamedTuple`: Table of specifications for a given partition k.
- `p::NamedTuple`: collection of key-value pairs where each key corresponds to a model parameter.
- `settings::Dict`: dictionary of settings containing configuration for the likelihood calculation.
- `fit_range`: array of arrays, defining the allowed energy ranges; e.g. `fit_range= [[1930,1950], [1970,1990], [2000,2050]]`.
- `bkg_units::String`: Specifies the units for the background index; available options are `"ckky"` (=counts/keV/kg/yr) or `"cFty"` (=counts/FWHM/t/yr).
"""
function build_likelihood_zero_obs_evts(
    part_k::NamedTuple,
    p::NamedTuple,
    settings::Dict,
    fit_range,
    bkg_units::String,
)

    ll_value = 0
    model_s_k, model_b_k = get_mu_s_b(p, part_k, 0, settings, fit_range, bkg_units)
    model_tot_k = model_b_k + model_s_k

    ll_value += -(model_tot_k + eps(model_tot_k))

    return ll_value
end


"""
    build_likelihood_per_partition(idx_part_with_events::Int,part_k::NamedTuple, events_k::Vector{Union{Float64}},p::NamedTuple,settings::Dict,bkg_shape::Symbol,fit_range,bkg_units::String)

Function which computes the likelihood for a single data partition k.

### Arguments
- `idx_part_with_events::Int`: index of the partition with the event.
- `part_k::NamedTuple`: Table of specifications for a given partition k.
- `events_k::Vector{Union{Float64}}`: vecotr of events (=energies) in the partition k.
- `p::NamedTuple`: collection of key-value pairs where each key corresponds to a model parameter.
- `settings::Dict`: dictionary of settings containing configuration for the likelihood calculation.
- `bkg_shape::Symbol`: Specifies the background shape; default is `:uniform`.
- `fit_range`: array of arrays, defining the allowed energy ranges; e.g. `fit_range= [[1930,1950], [1970,1990], [2000,2050]]`.
- `bkg_units::String`: Specifies the units for the background index; available options are `"ckky"` (=counts/keV/kg/yr) or `"cFty"` (=counts/FWHM/t/yr).
"""
function build_likelihood_per_partition(
    idx_part_with_events::Int,
    part_k::NamedTuple,
    events_k::Vector{Union{Float64}},
    p::NamedTuple,
    settings::Dict,
    bkg_shape::Symbol,
    fit_range,
    bkg_units::String,
)
    Qbb = ZeroNuFit.Constants.Qbb

    ll_value = 0

    model_s_k, model_b_k =
        get_mu_s_b(p, part_k, idx_part_with_events, settings, fit_range, bkg_units)

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
            reso, bias = ZeroNuFit.Utils.get_energy_scale_pars(
                part_k,
                p,
                settings,
                idx_part_with_events,
            )
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
    build_likelihood_looping_partitions(partitions::TypedTables.Table,events::Array{Vector{Float64}},part_event_index::Vector{Int},settings::Dict,sqrt_prior::Bool,s_max::Union{Float64,Nothing},fit_ranges,bkg_units::String;bkg_shape::Symbol=:uniformm)

Function to build the likelihood (a `DensityInterface.logfuncdensity` object) for the fit by looping over partitions.

### Arguments
- `partitions::TypedTables.Table`: table of partitions.
- `events::Array{Vector{Float64}}`: list of events (=energies) in each partition.
- `part_event_index::Vector{Int}`: index mapping events to the partitions.
- `settings::Dict`: A dictionary of settings containing configuration for the likelihood calculation.
- `sqrt_prior::Bool`: Whether to include the square root prior in the likelihood calculation. If `False`, a uniform prior is used.
- `s_max::Union{Float64, Nothing}`: A maximum value used for scaling the square root prior. If `Nothing`, no prior is applied.
- `fit_ranges`: The fitting ranges corresponding to the partitions.
- `bkg_shape::Symbol`: Specifies the background shape; default is `:uniform`.
- `bkg_units::String`: Specifies the units for the background index; available options are `"ckky"` (=counts/keV/kg/yr) or `"cFty"` (=counts/FWHM/t/yr).
"""
function build_likelihood_looping_partitions(
    partitions::TypedTables.Table,
    events::Array{Vector{Float64}},
    part_event_index::Vector{Int},
    settings::Dict,
    sqrt_prior::Bool,
    s_max::Union{Float64,Nothing},
    fit_ranges,
    bkg_units::String;
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
                        part_event_index[idx_k],
                        part_k,
                        events[idx_k],
                        p,
                        settings,
                        bkg_shape,
                        fit_ranges[part_k.fit_group],
                        bkg_units,
                    )
                else
                    # no events are there for a given partition
                    total_ll += build_likelihood_zero_obs_evts(
                        part_k,
                        p,
                        settings,
                        fit_ranges[part_k.fit_group],
                        bkg_units,
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
    generate_data(samples::BAT.DensitySampleVector,partitions::TypedTables.Table,part_event_index::Vector{Int},settings::Dict,fit_ranges,bkg_units::String;best_fit::Bool=false,seed=nothing,bkg_only=false)

Generates data from a posterior distribution.
This is based on the posterior predictive distributions. 
Given a model with some parameters `theta_i`, the posterior predictive distribution, or the distribution of data generated according to the posterior distribution of theta and the likelihood is:

```math
p(y|D) =int p(y|theta)p(theta|D)dtheta
```

Or in terms of sampling we first draw samples of `theta` from the posterior and then generate, datasets based on the likelihood.
We also give the options to fix the posterior distribution to the best fit, which is equivalent to the standard sampling methods.

### Arguments
- `samples::DensitySamplesVector`: samples of a past fit or a NamedTuple of best fit.
- `partitions::Table`: Table of the partition info.
- `part_event_index::Vector{Int}`: index for the parameters for partitions with events.
- `settings::Dict`: A dictionary of settings containing configuration for the likelihood calculation.
- `fit_ranges`: dictionary of energy ranges considered for the analysis.
- `best_fit::Bool`: True if you want to fix the paramaters to the best fit.
- `seed::Int`: random seed.
- `bkg_only::Bool`: True if we are using a model with background only.
- `bkg_units::String`: Specifies the units for the background index; available options are `"ckky"` (=counts/keV/kg/yr) or `"cFty"` (=counts/FWHM/t/yr).
"""
function generate_data(
    samples,#::BAT.DensitySampleVector,
    partitions::TypedTables.Table,
    part_event_index::Vector{Int},
    settings::Dict,
    fit_ranges,
    bkg_units::String;
    best_fit::Bool = false,
    seed = nothing,
    bkg_only = false,
)
    Qbb = ZeroNuFit.Constants.Qbb

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
            bkg_units,
        )

        n_s = rand(Poisson(model_s_k))
        n_b = rand(Poisson(model_b_k))
        events = ZeroNuFit.Utils.generate_disjoint_uniform_samples(
            n_b,
            fit_ranges[part_k.fit_group],
        )
        if (bkg_only == false)
            for i = 1:n_s

                reso, bias = ZeroNuFit.Utils.get_energy_scale_pars(
                    part_k,
                    p,
                    settings,
                    idx_part_with_events,
                )

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

Defines specific priors for signal and background contributions.

### Arguments
- `config`: input dictionary.
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

Builds the priors for use in the fit.

### Arguments
- `partitions::Table`: Table of the partition info.
- `part_event_index::Vector{Int}`: index for the parameters for partitions with events.
- `config`: input dictionary.
- `settings::Dict`: A dictionary of settings containing configuration for the likelihood calculation.
- `hierachical`: True if we use a hierarchical model.
- `hierachical_mode`: average mode of the hierarchical model.
- `hierachical_range`: range of the reference distribution used to correlated the background indexes in a hierarchical model.
- `bkg_shape::Symbol`: Specifies the background shape; default is `:uniform`.
- `shape_pars`: background shape parameters, e.g. slope for linear or exponential background modeling.
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
        :ε => [],
        :ω => [],
        :𝛥 => [],
    )

    for key in keys(distrB_multi)
        if config["bkg"]["units"] == "ckky"
            pretty_names[key] = string(key) * " [cts/keV/kg/yr]"
        end
        if config["bkg"]["units"] == "cFty"
            pretty_names[key] = string(key) * " [cts/FWHM/t/yr]"
        end
    end


    # create priors one by one
    ### SIGNAL PRIOR

    priors = OrderedDict()
    if (settings[:bkg_only] == false)
        priors[:S] = distrS
        @info "entered to add S prior"
    end

    # dictionary with info on the prior parameters
    nuisance_info =
        OrderedDict("α" => [], "αr" => [], "αb" => [], "ε" => [], "ω" => [], "𝛥" => [])

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
    if (settings[:eff_fixed] == false && settings[:eff_correlated] == false)
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

    if (settings[:energy_bias_fixed] == false && settings[:energy_bias_correlated] == true)

        list_names = partitions.energy_bias_name
        unique_list = unique(list_names)
        for name in unique_list
            priors[Symbol(name)] = Truncated(Normal(0, 1), -Inf, Inf)
            pretty_names[Symbol(name)] = L"\alpha_{b} (" * split(String(name), "_")[2] * ")"
            nuisance_info[string(name)] = [["combined", "", "", 0, 1, -Inf, Inf]]

        end
    end
    if (settings[:energy_bias_fixed] == false && settings[:energy_bias_correlated] == false)
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
                    # let's define some intervals in +-5σ 
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
    if (settings[:energy_res_fixed] == false && settings[:energy_res_correlated] == true)
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
    if (settings[:energy_res_fixed] == false && settings[:energy_res_correlated] == false)
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
end
