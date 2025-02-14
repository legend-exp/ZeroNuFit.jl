module Analysis

using Random, LinearAlgebra, Statistics, Distributions, StatsBase
using BAT, DensityInterface, IntervalSets
using TypedTables
using Plots
using Cuba
using SpecialFunctions

export run_analysis, retrieve_real_fit_results

using ZeroNuFit

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
            ZeroNuFit.Utils.save_generated_samples(samples, output_path)
            @info "...done!"
        end
    end

    @info "... now we save other useful results + config entries"
    ZeroNuFit.Utils.save_results_into_json(
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
        ZeroNuFit.Plotting.plot_correlation_matrix(
            samples,
            output_path,
            par_names = par_names,
            toy_idx = toy_idx,
        )
    end

    if config["light_output"] == false
        @info "... now we plot marginalized posteriors (and priors)"
        ZeroNuFit.Plotting.plot_marginal_distr(
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
    end

    if config["light_output"] == false
        @info "... now we plot 2D posterior"
        ZeroNuFit.Plotting.plot_two_dim_posteriors(
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
        ZeroNuFit.Plotting.plot_fit_and_data(
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

# function to run the unbinned fit
function run_analysis(config::Dict{String,Any}; output_path::String, toy_idx = nothing)
    """
    Function which handeles running analysis
    Parameters:
    ----------
        config::Dict{String,Any} the fit configuration
        output_path::String (keyword) the path to the output files folder
    """
    @info "You entered into src/ZeroNuFit.jl"

    part_event_index, events, partitions, fit_ranges =
        ZeroNuFit.Utils.get_partitions_events(config)
    # check if you want to overwrite the fit; if no results are present, then fit data
    if config["overwrite"] == true ||
       !isfile(joinpath(config["output_path"], "mcmc_files/samples.h5"))
        @info "... now we run a fit"

        if config["overwrite"] == true
            @info "OVERWRITING THE PREVIOUS FIT!"
        end

        samples, prior, par_names = ZeroNuFit.Likelihood.run_fit_over_partitions(
            partitions,
            events,
            part_event_index,
            config,
            fit_ranges,
        )
        @info "fit ran succesfully"
    else
        @info "... we load already existing fit results"
        samples = bat_read(joinpath(config["output_path"], "mcmc_files/samples.h5")).result
        prior, _, _, par_names, nuisance_info = ZeroNuFit.Likelihood.get_stat_blocks(
            partitions,
            events,
            part_event_index,
            fit_ranges,
            config = config,
            bkg_only = config["bkg_only"],
        )
    end

    # let's save
    @info samples
    @info bat_report(samples)
    _, _, posterior, _, nuisance_info = ZeroNuFit.Likelihood.get_stat_blocks(
        partitions,
        events,
        part_event_index,
        fit_ranges,
        config = config,
        bkg_only = config["bkg_only"],
    )
    save_outputs(
        partitions,
        events,
        part_event_index,
        samples,
        posterior,
        nuisance_info,
        config,
        output_path,
        fit_ranges,
        priors = prior,
        par_names = par_names,
        toy_idx = toy_idx,
    )

    return

end


function retrieve_real_fit_results(config::Dict{String,Any})
    """
    Function which handeles generating of fake data
    Parameters:
    ----------
        config::Dict{String,Any} the fit configuration
        output_path::String (keyword) the path to the output files folder
    """

    @info "You entered into src/ZeroNuFit.jl"

    @info"Let's retrieve some partitions ..."
    partitions = nothing
    first = true
    @info config["partitions"]
    for part_path in config["partitions"]

        part_tmp, fit_groups = ZeroNuFit.Utils.get_partitions_new(part_path)
        if (first)
            partitions = part_tmp
            first = false
        else
            partitions = vcat(partitions, part_tmp)
        end
    end
    display(partitions)
    @info "... load events"
    events_multi = []
    for event_path in config["events"]
        append!(events_multi, [ZeroNuFit.Utils.get_events(event_path, partitions)])
    end

    events = Array{Vector{Float64}}(undef, length(partitions))
    for i = 1:length(partitions)

        arr_tmp = Vector{Float64}()
        for sub in events_multi
            if (sub[i] != Float64[])

                append!(arr_tmp, sub[i])
            end
        end

        events[i] = arr_tmp
    end
    @debug events

    @info "get which partitions have events"
    part_event_index = ZeroNuFit.Utils.get_partition_event_index(events, partitions)

    # let's retrieve the old data
    samples = bat_read(joinpath(config["output_path"], "mcmc_files/samples.h5")).result

    return samples, partitions, part_event_index

end
end
