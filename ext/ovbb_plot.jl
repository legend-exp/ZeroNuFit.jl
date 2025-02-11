using Pkg
Pkg.activate(".")
using StatsBase: weights, Histogram
using PyCall
@pyimport matplotlib.pyplot as plt
using BAT
import HDF5
using IntervalSets: (..)

include("../src/utils.jl")
include("../src/fitting.jl")
include("../src/constants.jl")
include("../src/likelihood.jl")

# this function calculates the observable counts density at a certain energy
# and given the samples. "density" means that there is no integral over energy
# (i.e. no multiplication by ΔE)
function counts_density_obs(partitions, energy, S, α, pars, fit_ranges; shape=nothing, bkg = true, signal = true)
    μk_tot = 0

    for P in partitions
        ϵk = P.eff_tot
        ϵk_sigma = P.eff_tot_sigma
        μsk = S * log(2) * constants.N_A * P.exposure .* (ϵk .+ α * ϵk_sigma) / constants.m_76

        if shape == "linear"
            b_name = P.bkg_name
            slope_name = Symbol(String(P.bkg_name) * "_slope")
            slope = pars[slope_name]

            fit_range = fit_ranges[P.fit_group]
            range_l, range_h = get_range(fit_range)
            center = range_l[1]
            delta = range_h[end] - range_l[1]

            μbk = pars[b_name] * P.exposure * (1 + slope * (energy - center) / delta) 

        elseif shape == "exponential"
            b_name = P.bkg_name
            slope_name = Symbol(String(P.bkg_name) * "_slope")
            slope = pars[slope_name]

            fit_range = fit_ranges[P.fit_group]
            range_l, range_h = get_range(fit_range)
            center = range_l[1]
            R = p[Symbol(string(b_name) * "_slope")]
            delta = range_h[end] - range_l[1]
            Rt = R / delta
            μbk = pars[b_name] * P.exposure * exp_stable((energy - center) * Rt)
            
        else
            b_name = P.bkg_name
            μbk = pars[b_name] * P.exposure
        end

        Δk = P.bias
        Δk_sigma = P.bias_sigma
        σk = P.width
        σk_sigma = P.width_sigma

        # NOTE: this assumes that the posterior for bias/resolution in "get_signal_pdf" is equal to the prior
        μk =
            (bkg ? μbk : 0) .+ (
                signal ?
                μsk .* get_signal_pdf(energy, constants.Qbb, P) * constants.sig_units : 0
            )
        μk_tot += μk
    end

    return μk_tot
end

# this function calculates the posterior of the counts density/(kg yr) at a certain energy
function counts_density_kgyr_l200_post(energy, samples, partitions, fit_ranges; kwargs...)
    @info "calculating posterior for energy $energy"
    f = s -> counts_density_obs(partitions, energy, s.S, s.αe_all, s, fit_ranges; kwargs...)
    return broadcast(f, samples.v) ./ sum(partitions.exposure)
end

# this function calculates the posterior of the counts density/(kg yr) at a certain energy, for bkg fixed at the best fit value
function counts_density_kgyr_l200_bkg_mode(energy, mode, samples, partitions, fit_ranges; kwargs...)
    @info "calculating posterior for energy $energy"
    f = s -> counts_density_obs(partitions, energy, s.S, s.αe_all, mode, fit_ranges; kwargs...)
    return broadcast(f, samples.v) ./ sum(partitions.exposure)
end

# this function calculates the counts density/(kg yr) at the best fit value
function counts_density_kgyr_l200_mode(energy, best_fit, partitions, fit_ranges; kwargs...)
    n_cts = counts_density_obs(
        partitions,
        energy,
        best_fit.S,
        best_fit.αe_all,
        best_fit, 
        fit_ranges;
        kwargs...,
    )
    return n_cts ./ sum(partitions.exposure)
end


function plot_l200_result(samples, part_event_index, events, partitions, fit_ranges, config)

    # we don't need all samples
    _samples = samples[1:20_000]
    _weights = weights(_samples.weight)

    # for the background, we just need one (random) energy point to compute the 68% smallest CI if linear;
    # otherwise, we have to compute it for at least 2 values for the linear case;
    # or multiple energy values for the exponential case

    # let's define a flag for the bkg shape
    if "shape" in keys(config["bkg"])
        bkg_shape = config["bkg"]["shape"]["name"]
    else
        bkg_shape = nothing
    end

    # ... 2 points are enough for the linear case (start & end of the analysis window)
    energies_all_window = 1930:20:2190  #[1930, 1960, 2000, 2050, 2190]
    # ...more points are needed for the exponential case
    if bkg_shape == "exponential"
        energies_all_window = 1930:100:2190 # INCREASE!!!
    end
    # if flat, evaluate 68% band at one energy value (eg 2000 keV)
    if bkg_shape == nothing
        posterior = counts_density_kgyr_l200_post(2000, _samples, partitions, fit_ranges, shape=bkg_shape, signal = false)
        _int = BAT.smallest_credible_intervals(posterior, _weights, nsigma_equivalent = 1)
        length(_int) != 1 && @warn "[2000 keV] 68% interval is disjoint", _int
        b_68 = first(_int)
    # if not flat, ...
    else
        # now we have a vector with as many entries as the energies we selected
        b_68 = Vector(undef, length(energies_all_window))
        # let's speed this up with multi-threading
        Threads.@threads for i = 1:length(energies_all_window)
            E = energies_all_window[i]
            posterior = counts_density_kgyr_l200_post(E, _samples, partitions, fit_ranges, shape=bkg_shape, signal = false)
            _int = BAT.smallest_credible_intervals(posterior, _weights, nsigma_equivalent = 1)
            length(_int) != 1 && @warn "[$E keV] 68% interval is disjoint: ", _int
            b_68[i] = first(_int)
        end
    end   
    
    # best fit results    
    _, _, posterior, _, _ = get_stat_blocks(
        partitions,
        events,
        part_event_index,
        fit_ranges,
        config = config,
        bkg_only = config["bkg_only"],
    )
    _, best_fit_pars = get_global_mode(samples, posterior)

    # let's use a narrow region around Qbb for the 90% CI on the signal
    energies = 2034.0:0.2:2044.0
    s_90 = Vector(undef, length(energies))

    b_mode_all_window = Vector(undef, length(energies_all_window))
    
    if bkg_shape == nothing
        b_mode = counts_density_kgyr_l200_mode(2000, best_fit_pars, partitions, fit_ranges, shape=bkg_shape, signal = false)
    else
        # we use the same energies of the signal (this will be used when plotting the signal)
        b_mode = Vector(undef, length(energies))
        Threads.@threads for i = 1:length(energies)
            E = energies[i]
            b_mode[i] = counts_density_kgyr_l200_mode(E, best_fit_pars, partitions, fit_ranges, shape=bkg_shape, signal = false)
        end

        # we use all energies for the mode all over the fit window
        Threads.@threads for i = 1:length(energies_all_window)
            E = energies_all_window[i]
            b_mode_all_window[i] = counts_density_kgyr_l200_mode(E, best_fit_pars, partitions, fit_ranges, shape=bkg_shape, signal = false)
        end
    end

    # extract signal upper limits
    Threads.@threads for i = 1:length(energies)
        E = energies[i]
        posterior = counts_density_kgyr_l200_bkg_mode(E, best_fit_pars, _samples, partitions, fit_ranges, shape=bkg_shape, bkg=false)
        if bkg_shape == nothing
            s_90[i] = b_mode .. quantile(posterior, _weights, 0.9)
        else
            s_90[i] = b_mode[i] .. quantile(posterior, _weights, 0.9)
        end
    end

    fig, ax = plt.subplots(figsize = (5, 2.5))

    # ...background 68%
    if bkg_shape == nothing
        band = plt.fill_between(
            [1930, 2190],
            fill(b_68.left, 2),
            fill(b_68.right, 2),
            linewidth = 0,
            alpha = 0.5,
            color = "#228833",
            label = "68% C.I.",
            zorder = 0,
        )
    else
        left_values = [i.left for i in b_68]
        right_values = [i.right for i in b_68]
        band = plt.fill_between(energies_all_window, left_values, right_values, color="#228833", label = "68% C.I.", alpha=0.5, linewidth=0, zorder=0)
    end
    
    # ...background best fit
    if bkg_shape == nothing
        line, = plt.plot([1930, 2190], [b_mode, b_mode], color = "#228833", zorder = 1)
    else
        line, = plt.plot(energies_all_window, b_mode_all_window, color = "#228833", zorder = 1)
    end

    # ...signal 90%
    limit_area = plt.fill_between(
        energies,
        getproperty.(s_90, :left),
        getproperty.(s_90, :right),
        linewidth = 0,
        color = "#4477AA",
        label = "90% C.I.",
    )

    # plot the data
    unbinned = true

    # get energies of surviving events
    unflat_x = vcat([y for y in events if !isempty(y)])
    x = vcat(unflat_x...)
    data_art = nothing
    exp_tot = sum(partitions.exposure)

    if unbinned
        data_art, stemlines, baseline =
            ax.stem(x, fill(1e-3, length(x)), "#CC3311", basefmt = "none")
        data_art.set_markersize(1)
        # NOTE: comment the following if the stem has a length!
        data_art.set_zorder(99)
        data_art.set_clip_on(false)

        stemlines.set_linewidth(0.6)
    else
        w = 1 / exp_tot
        _, _, data_art =
            plt.hist(x, weights = fill(w, length(x)), bins = 1930:1:2190, color = "#CC3311")
    end

    ax.set_ylim(1E-4, 2)
    ax.set_yscale("log")
    plt.axvspan(2099, 2109, facecolor="black", alpha=0.4)
    plt.axvspan(2114, 2124, facecolor="black", alpha=0.4)
    ax.legend(
        [data_art, (band, line), limit_area],
        [
            string(raw"LEGEND-200 data [", round(exp_tot, digits = 1), " kg yr]"),
            raw"Background best fit and 68% C.I. interval",
            string(raw"$T^{0\nu}_{1/2} > $", limit_half_life, raw"$\times 10^{26}$ yr [90% C.I.]"),
        ],
        loc = "upper left",
        fontsize = "small",
    )

    output = "../ZeroNuFit-dev/output_correct_inputs_uncorrelated_nuisances/v3/fit_9_l200_uniform_1BI_CorrEff"#config["output"]
    fit_results = JSON.parsefile(joinpath(output, "mcmc_files/fit_results.json"))["quantile90"]
    limit_signal = fit_results["S"]
    limit_half_life = round(10/(limit_signal), digits=1)
    @info "limit_signal -> ", limit_signal
    @info "limit_half_life -> ", limit_half_life
    ax.set_xlim(1930, 2190)
    ax.set_xlabel("Energy [keV]", loc = "right")
    ax.set_ylabel("Counts / (keV kg yr)")
    plt.tight_layout()
    plt.savefig("l200_result.pdf", bbox_inches = "tight")
end


function get_total_bi(samples, part_event_index, events, partitions, fit_ranges, config)
    energy = 2000

    function get_bi(partitions, energy, pars)
        μb_tot = 0        
        for P in partitions
            b_name = P.bkg_name
            μbk = pars[b_name] * P.exposure
            μb_tot+= μbk
        end
        return μb_tot
    end

    f = samples -> get_bi(partitions, energy, samples)
    bi_tot = broadcast(f, samples.v)  ./ sum(partitions.exposure)
    h5write("bi_tot_all.h5", "background", bi_tot)
end


samples = bat_read(joinpath("../ZeroNuFit-dev/output_correct_inputs_uncorrelated_nuisances/v3/fit_9_l200_uniform_1BI_CorrEff_linear", "mcmc_files/samples.h5")).result
config_file = "legend-0vbb-config/config/v3/fit_9_l200_uniform_1BI_CorrEff_linear.json"
config = JSON.parsefile(config_file)
part_event_index, events, partitions, fit_ranges = get_partitions_events(config)
plot_l200_result(samples, part_event_index, events, partitions, fit_ranges, config)
#get_total_bi(samples, part_event_index, events, partitions, fit_ranges, config)