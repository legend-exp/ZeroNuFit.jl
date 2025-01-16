using PropertyFunctions
using JSON
using Logging
using Random, LinearAlgebra, Statistics, Distributions, StatsBase
using PropDicts
using FilePathsBase
using DataStructures
using PropDicts
using Tables
using TypedTables
using Optim
using FileIO
import JLD2
import HDF5

function check_key(config, k)
    if !(k in keys(config))
        @error "'$k' not in config, exit here"
        exit(-1)
    end
end

function get_settings(config)

    check_key(config, "nuisance")
    check_key(config["nuisance"], "energy_scale")
    check_key(config["nuisance"]["energy_scale"], "fixed")

    settings = Dict()
    settings[:energy_scale_fixed] = config["nuisance"]["energy_scale"]["fixed"]
    settings[:energy_scale_correlated] = config["nuisance"]["energy_scale"]["correlated"]
    settings[:eff_fixed] = config["nuisance"]["efficiency"]["fixed"]
    settings[:eff_correlated] = config["nuisance"]["efficiency"]["correlated"]
    settings[:bkg_only] = config["bkg_only"]

    return settings
end


"""
    get_partitions_new(part_path::String)

Get the partition info from a JSON file and save to a Table
"""
function get_partitions_new(part_path::String)
    if !isfile(part_path)
        @error "Error: file $part_path does not exist!"
        exit(-1)
    end
    part_data_json = JSON.parsefile(part_path, dicttype = DataStructures.OrderedDict)

    fit_groups = part_data_json["fit_groups"]

    list_groups = collect(keys(fit_groups))
    k = keys(part_data_json["partitions"][list_groups[1]][1])
    arrays = Dict()
    for key in k
        arrays[key] = []

    end
    arrays["fit_group"] = []
    arrays["signal_par_name"] = []
    arrays["bkg_par_name"] = []
    arrays["eff_par_name"] = []
    arrays["energy_reso_name"] = []
    arrays["energy_bias_name"] = []
    arrays["frac"] = []
    arrays["tau"] = []
    arrays["sigma"] = []

    fit_ranges = OrderedDict()
    for fit_group in keys(part_data_json["partitions"])

        fit_ranges[fit_group] = part_data_json["fit_groups"][fit_group]["range"]
        for part in part_data_json["partitions"][fit_group]
            for key in k
                if key in ["frac", "tau", "sigma"]
                    continue # separate treatment
                end
                append!(arrays[key], [part[key]])
            end
            for key in ["frac", "tau", "sigma"]
                if key in k
                    append!(arrays[key], [part[key]])
                else
                    append!(arrays[key], [nothing])
                end
            end
            append!(arrays["fit_group"], [fit_group])
            append!(
                arrays["bkg_par_name"],
                [Symbol(part_data_json["fit_groups"][fit_group]["bkg_name"])],
            )

            ## defaults to 'gaussian'
            if haskey(part_data_json["fit_groups"][fit_group], "signal_name")
                append!(
                    arrays["signal_par_name"],
                    [Symbol(part_data_json["fit_groups"][fit_group]["signal_name"])],
                )
            else
                append!(arrays["signal_par_name"], [Symbol("gaussian")])
            end

            ## defaults to 'all'
            if (haskey("efficiency_group_name", part_data_json["fit_groups"][fit_group]))
                append!(
                    arrays["eff_par_name"],
                    [
                        "αe_" * Symbol(
                            part_data_json["fit_groups"][fit_group]["efficiency_group_name"],
                        ),
                    ],
                )
            else
                append!(arrays["eff_par_name"], [:αe_all])
            end

            ## defaults to 'all'
            if (haskey("energy_scale_group_name", part_data_json["fit_groups"][fit_group]))
                append!(
                    arrays["energy_reso_name"],
                    [
                        Symbol(
                            "αr_" *
                            part_data_json["fit_groups"][fit_group]["energy_scale_group_name"],
                        ),
                    ],
                )
                append!(
                    arrays["energy_bias_name"],
                    [
                        Symbol(
                            "αb_" *
                            part_data_json["fit_groups"][fit_group]["energy_scale_group_name"],
                        ),
                    ],
                )

            else
                append!(arrays["energy_reso_name"], [:αr_all])
                append!(arrays["energy_bias_name"], [:αb_all])

            end

        end

    end

    #TODO: find a way to make this not hardcoded
    tab = Table(
        experiment = Array(arrays["experiment"]),
        fit_group = Array(arrays["fit_group"]),
        bkg_name = Array(arrays["bkg_par_name"]),
        signal_name = Array(arrays["signal_par_name"]),
        energy_reso_name = Array(arrays["energy_reso_name"]),
        energy_bias_name = Array(arrays["energy_bias_name"]),
        eff_name = Array(arrays["eff_par_name"]),
        detector = Array(arrays["detector"]),
        part_name = Array(arrays["part_name"]),
        start_ts = Array(arrays["start_ts"]),
        end_ts = Array(arrays["end_ts"]),
        eff_tot = Array(arrays["eff_tot"]),
        eff_tot_sigma = Array(arrays["eff_tot_sigma"]),
        width = Array(arrays["width"]),
        width_sigma = Array(arrays["width_sigma"]),
        exposure = Array(arrays["exposure"]),
        bias = Array(arrays["bias"]),
        bias_sigma = Array(arrays["bias_sigma"]),
        frac = Array(arrays["frac"]),
        tau = Array(arrays["tau"]),
        sigma = Array(arrays["sigma"]),
    )
    return tab, fit_groups, fit_ranges
end


"""
    event_is_contained(event::Float64, fit_ranges)
    
Returns true if the event is contained at least in one of the energy ranges
"""
function event_is_contained(event::Float64, fit_range)::Bool
    flag = false
    for range_pair in fit_range
        if range_pair[1] <= event <= range_pair[2]
            flag = true
        end
    end
    return flag
end


"""
    get_partitions_events(config::Dict{String, Any})
    
Get partition, event, and fit range info from input JSON files
"""
function get_partitions_events(config::Dict{String,Any})

    @info"... retrieve some partitions"
    partitions, fit_ranges = get_partitions(config)
    display(partitions)
    @info "... load events"
    events_multi = []
    for event_path in config["events"]
        append!(events_multi, [get_events(event_path, partitions)])
    end

    events = Array{Vector{Float64}}(undef, length(partitions))
    for i = 1:length(partitions)

        # get fit group (each one has its own fit range) 
        fit_group = partitions[i].fit_group
        fit_range = fit_ranges[fit_group]
        arr_tmp = Vector{Float64}()
        for sub in events_multi
            if (sub[i] != Float64[])
                # check if the event is contained in the given fit range
                if event_is_contained(sub[i][1], fit_range)
                    append!(arr_tmp, sub[i])
                end
            end
        end

        events[i] = arr_tmp
    end
    @debug events

    @info "... get which partitions have events"
    part_event_index = get_partition_event_index(events, partitions)

    return part_event_index, events, partitions, fit_ranges

end


"""
    get_partition_event_index(events::Array{Vector{Float64}},partitions::TypedTables.Table)::Vector{Int}

Gets an object descirbing if a partition has an event and giving them indexes
This creates a vector where

- V[i]=0 if partition i has no events
- V[i]=idx if partition i has events

where the index counts the number of partitions with index<=i with, 
events and corresponds to the index of the parameters.
"""
function get_partition_event_index(
    events::Array{Vector{Float64}},
    partitions::TypedTables.Table,
)::Vector{Int}
    output = Vector{Int}(undef, length(partitions))
    counter = 1
    for (idx, part) in enumerate(partitions)
        if (events[idx] != Any[])
            output[idx] = counter
            counter += 1
        else
            output[idx] = 0
        end
    end

    return output
end


"""
    get_partitions(config::Dict{String, Any})

Get the partition info from a JSON file 
"""
function get_partitions(config::Dict{String,Any})

    partitions = nothing
    first = true
    fit_ranges = nothing

    check_key(config, "partitions")
    check_key(config, "events")

    for part_path in config["partitions"]

        part_tmp, fit_groups, fit_range = get_partitions_new(part_path)
        if (first)
            partitions = part_tmp
            first = false
            fit_ranges = fit_range
        else
            partitions = vcat(partitions, part_tmp)
            merge!(fit_ranges, fit_range)
        end
    end
    return partitions, fit_ranges
end


"""
    get_events(event_path::String,partitions)::Array{Vector{Float64}}

Get the event info from a JSON file and save to a Table
"""
function get_events(event_path::String, partitions)::Array{Vector{Float64}}
    if !isfile(event_path)
        @error "Error: file $event_path does not exist!"
        exit(-1)
    end
    event_json = JSON.parsefile(event_path, dicttype = DataStructures.OrderedDict)
    events = Array{Vector{Float64}}(undef, length(partitions))
    for (idx, part) in enumerate(partitions)
        events[idx] = Vector{Float64}[]
    end

    for event in event_json["events"]
        found = false

        for (idx, part) in enumerate(partitions)
            if (
                part.experiment == event["experiment"] &&
                part.detector == event["detector"] &&
                event["timestamp"] < part.end_ts + 1 &&
                event["timestamp"] > part.start_ts - 1
            )
                append!(events[idx], Vector{Float64}([Float64(event["energy"])]))
                found = true
            end
        end

        if (found == false)
            @error event "has no partition"
            exit(-1)
        end
    end

    return events

end

"""
    get_efficiency(p::NamedTuple,part_k::NamedTuple,idx_part_with_events::Int,settings::Dict)

Get the efficiency 
"""
function get_efficiency(
    p::NamedTuple,
    part_k::NamedTuple,
    idx_part_with_events::Int,
    settings::Dict,
)
    eff = nothing

    # if background only fit, then there is no need to retrieve any efficiency 
    # (enters only in the signal-related term) 
    if settings[:bkg_only] == true
        return 0
    end

    # CORRELATED efficiency (eff is described by one global alpha parameter)
    if (settings[:eff_correlated] == true)
        eff_group = part_k.eff_name
        eff = part_k.eff_tot + p[eff_group] * part_k.eff_tot_sigma

        # UNCORRELATED efficiency (eff follows a pdf, different for each partition)
    elseif (
        idx_part_with_events != 0 &&
        settings[:eff_correlated] == false &&
        settings[:eff_fixed] == false
    )
        eff = p.ε[idx_part_with_events]

        # FIXED efficiency
    else
        eff = part_k.eff_tot
    end

    return eff
end


""" 
    get_energy_scale_pars(part_k::NamedTuple,p::NamedTuple,settings::Dict,idx_part_with_events)

Get the resolution and bias
"""
function get_energy_scale_pars(
    part_k::NamedTuple,
    p::NamedTuple,
    settings::Dict,
    idx_part_with_events,
)
    # either FIXED energy pars OR 0 events in the partition
    if (settings[:energy_scale_fixed] == true || idx_part_with_events == 0)
        reso = part_k.width
        bias = part_k.bias
        # CORRELATED energy pars (reso and bias are described by one global alpha parameter each)
    elseif (settings[:energy_scale_correlated] == true)
        energy_reso_group = part_k.energy_reso_name
        energy_bias_group = part_k.energy_bias_name
        reso = part_k.width + p[energy_reso_group] * part_k.width_sigma
        bias = part_k.bias + p[energy_bias_group] * part_k.bias_sigma
        # UNCORRELATED energy pars (reso and bias follow a pdf each, different for each partition)
    else
        reso = p.ω[idx_part_with_events]
        bias = p.𝛥[idx_part_with_events]
    end

    # convert into Float
    return reso * 1.0, bias * 1.0

end


## sampling 
function inverse_uniform_cdf(
    p,
    fit_range::Union{Vector{Vector{Int}},Vector{Vector{Float64}}},
)
    range_l, range_h = get_range(fit_range)
    delta = sum(range_h .- range_l)

    cumulative_prob = 0.0
    cumulative_range = 0.0

    # handle the p=1 case
    if p == 1
        res = range_h[end]
    else
        for j = 1:length(fit_range)
            interval_width = range_h[j] - range_l[j]
            interval_prob = interval_width / delta

            # if p is within the current interval range
            if cumulative_prob + interval_prob >= p
                # scale the probability to the current interval
                interval_p = (p - cumulative_prob) / interval_prob
                res = range_l[j] + interval_p * interval_width
                break
            end

            # accumulate the probability and range covered so far
            cumulative_prob += interval_prob
            cumulative_range += interval_width
        end
    end

    return res
end


function generate_disjoint_uniform_samples(n, fit_range; seed = nothing)
    # fix the seed (if provided)
    if seed !== nothing
        Random.seed!(seed)
    end
    rands = rand(n)
    res = [inverse_uniform_cdf(rand, fit_range) for rand in rands]
    return res
end


"""
    save_generated_samples(samples,output)

Function which saves sampling results
"""
function save_generated_samples(samples, output)
    FileIO.save(joinpath(output, "mcmc_files/samples.jld2"), Dict("samples" => samples))
    bat_write(joinpath(output, "mcmc_files/samples.h5"), samples)
end


"""
    get_global_mode(samples, posterior)

Function which retrieves global mode and a refined estimate of it
"""
function get_global_mode(samples, posterior)
    global_modes = BAT.mode(samples)
    # more refined estimate 
    findmode_result = bat_findmode(
        posterior,
        OptimAlg(optalg = Optim.NelderMead(), init = ExplicitInit([global_modes])),
    )
    return global_modes, findmode_result.result
end

"""
    get_marginalized_mode(samples, par)

Function which retrieves marginalized mode as the highest bin of the posterior, rebinned with 250 bins
"""
function get_marginalized_mode(samples, par)
    post = get_par_posterior(samples, par, idx = nothing)
    post_numeric = Float64.(post)
    xmin = minimum(post)
    xmax = maximum(post)
    nbin = 250
    delta = (xmax - xmin) / nbin
    hist = fit(Histogram, post_numeric, xmin:delta:xmax)
    max_bin_idx = argmax(hist.weights)
    mode_value = hist.edges[1][max_bin_idx]
    return mode_value
end


"""
    save_results_into_json(samples,posterior,nuisance_info,config,output;par_names=nothing,toy_idx=nothing)

Function which saves results from the fit and copies the input config (for any future need)
"""
function save_results_into_json(
    samples,
    posterior,
    nuisance_info,
    config,
    output;
    par_names = nothing,
    toy_idx = nothing,
)

    global_modes, refined_global_modes = get_global_mode(samples, posterior)

    # save partitions info for nuisance parameters
    nuisance_dict = Dict{String,Vector{Dict{String,Any}}}()
    for (key, values) in pairs(nuisance_info)
        # skip parameters that were not used in the fit
        if values == []
            continue
        end
        # initialize an empty array for each main key (eff, bias, res)
        nuisance_dict[key] = []

        for idx = 1:length(values)

            inner_dict = Dict(
                "experiment" => values[idx][1],
                "partition" => values[idx][2],
                "detector" => values[idx][3],
                "prior_mu" => values[idx][4],
                "prior_sigma" => values[idx][5],
                "prior_low_bound" => values[idx][6] == -Inf ? "-Inf" : values[idx][6],
                "prior_upp_bound" => values[idx][7] == Inf ? "Inf" : values[idx][7],
            )
            push!(nuisance_dict[key], inner_dict)
        end
    end

    # default marginalized modes
    ltmp = NullLogger()
    marginalized_modes = 0
    with_logger(ltmp) do
        marginalized_modes = BAT.bat_marginalmode(samples).result
    end

    # marginalized mode from binned histogram (250 bins)
    unshaped_samples, f_flatten = bat_transform(Vector, samples)
    first_sample = unshaped_samples.v[1]
    free_pars = keys(first_sample)
    marginalized_modes_highest_bin = Dict()
    ct = 1
    for (idx, par) in enumerate(keys(marginalized_modes))
        # parameters with 1 entry only
        if length(marginalized_modes[par]) == 1
            marginalized_modes_highest_bin[string(par)] =
                get_marginalized_mode(samples, par)
            ct += 1
            # parameters with more than 1 entry
        else
            marginalized_modes_highest_bin[string(par)] = []
            for entry in marginalized_modes[par]
                mode_value = get_marginalized_mode(unshaped_samples, free_pars[ct])
                ct += 1
                append!(marginalized_modes_highest_bin[string(par)], mode_value)
            end
        end
    end

    mean = BAT.mean(samples)
    stddev = BAT.std(samples)

    ci_68 = BAT.smallest_credible_intervals(samples, nsigma_equivalent = 1)
    ci_90 = BAT.smallest_credible_intervals(samples, nsigma_equivalent = 1.64)
    ci_95 = BAT.smallest_credible_intervals(samples, nsigma_equivalent = 2)
    ci_99 = BAT.smallest_credible_intervals(samples, nsigma_equivalent = 3)

    quantile90 = Statistics.quantile(samples, 0.9)

    data = Dict(
        "mean" => mean,
        "stddev" => stddev,
        "global_modes" => global_modes,
        "refined_global_modes" => refined_global_modes,
        "marginalized_modes" => marginalized_modes,
        "marginalized_modes_highest_bin" => marginalized_modes_highest_bin,
        "ci_68" => ci_68,
        "ci_90" => ci_90,
        "ci_95" => ci_95,
        "ci_99" => ci_99,
        "quantile90" => quantile90,
        "config" => config,
        "nuisance_info" => nuisance_dict,
    )

    json_string = JSON.json(data, 4)

    if toy_idx == nothing
        open(joinpath(output, "mcmc_files/fit_results.json"), "w") do file
            write(file, json_string)
        end
    else
        open(joinpath(output, "mcmc_files/fit_results_$(toy_idx).json"), "w") do file
            write(file, json_string)
        end
    end
end


"""
    save_outputs(partitions, events, part_event_index, samples, posterior, nuisance_info, config, output_path, fit_ranges;priors=nothing,par_names=nothing,toy_idx=nothing)

Function to plot and save results, as well as inputs
"""
function save_outputs(
    partitions,
    events,
    part_event_index,
    samples,
    posterior,
    nuisance_info,
    config,
    output_path,
    fit_ranges;
    priors = nothing,
    par_names = nothing,
    toy_idx = nothing,
)
    if (haskey(config["bkg"], "correlated")) &
       (config["bkg"]["correlated"]["mode"] != "none")
        hier = true
    else
        hier = false
    end
    if (config["signal"]["prior"] == "sqrt")
        sqrt_prior = true
        s_max = config["signal"]["upper_bound"]
    else
        sqrt_prior = false
        s_max = nothing
    end
    first_sample = samples.v[1]
    free_pars = keys(first_sample) # in format (:B, :S, ...) 
    @info "... these are the parameters that were included: ", free_pars

    @info "... now we save samples (untouched if we do not want to overwrite)"
    if config["light_output"] == false
        if config["overwrite"] == true ||
           !isfile(joinpath(config["output_path"], "mcmc_files/samples.h5"))
            save_generated_samples(samples, output_path)
            @info "...done!"
        end
    end

    @info "... now we save other useful results + config entries"
    save_results_into_json(
        samples,
        posterior,
        nuisance_info,
        config,
        output_path,
        par_names = par_names,
        toy_idx = toy_idx,
    )
    @info "...done!"

    if config["light_output"] == false
        plot_correlation_matrix(
            samples,
            output_path,
            par_names = par_names,
            toy_idx = toy_idx,
        )
    end

    @info "... now we plot marginalized posteriors (and priors)"
    plot_marginal_distr(
        partitions,
        samples,
        free_pars,
        output_path,
        priors = priors,
        par_names = par_names,
        plot_config = config["plot"],
        s_max = s_max,
        sqrt_prior = sqrt_prior,
        hier = hier,
        toy_idx = toy_idx,
    )

    if config["light_output"] == false
        @info "... now we plot 2D posterior"
        plot_two_dim_posteriors(
            samples,
            free_pars,
            output_path,
            par_names = par_names,
            toy_idx = toy_idx,
        )
        @info "...done!"
    end

    if config["plot"]["bandfit_and_data"] || config["plot"]["fit_and_data"]
        @info "... now we plot fit & data"
        plot_fit_and_data(
            partitions,
            events,
            part_event_index,
            samples,
            posterior,
            free_pars,
            output_path,
            config,
            fit_ranges,
            toy_idx = toy_idx,
        )
        @info "...done!"
    end

end
