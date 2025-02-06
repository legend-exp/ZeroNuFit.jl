using Pkg
Pkg.activate(".")
using StatsBase: weights, Histogram
using PyCall
@pyimport matplotlib.pyplot as plt
using BAT
import HDF5
using IntervalSets: (..)

include("utils.jl")
include("fitting.jl")
include("constants.jl")
include("likelihood.jl")

# this function calculates the observable counts density at a certain energy
# and given the samples. "density" means that there is no integral over energy
# (i.e. no multiplication by ΔE)
function counts_density_obs(partitions, energy, S, B, α; bkg=true, signal=true)
    μk_tot = 0

    for P in partitions
        ϵk = P.eff_tot
        ϵk_sigma = P.eff_tot_sigma

        μsk = S * log(2) * constants.N_A * P.exposure .* (ϵk .+ α * ϵk_sigma) / constants.m_76

        # TO-DO: this has to be generalized to the linear and exponential case (removing the normalization factor - refactor is needed)
        μbk = B * P.exposure

        Δk = P.bias
        Δk_sigma = P.bias_sigma
        σk = P.width
        σk_sigma = P.width_sigma

        # NOTE: this assumes that the posterior for bias/resolution in "get_signal_pdf" is equal to the prior
        μk = (bkg ? μbk : 0) .+ (signal ? μsk .* get_signal_pdf(energy, constants.Qbb, P) * constants.sig_units : 0)
        μk_tot += μk
    end

    return μk_tot
end

# this function calculates the posterior of the counts density/(kg yr) at a certain energy
function counts_density_kgyr_l200_post(energy, samples, partitions; kwargs...)
    @info "calculating posterior for energy $energy"
    f = s -> counts_density_obs(partitions, energy, s.S, s.B_l200a_all, s.αe_all; kwargs...)
    return broadcast(f, samples.v) ./ sum(partitions.exposure)
end

# this function calculates the posterior of the counts density/(kg yr) at a certain energy, for bkg fixed at the best fit value
function counts_density_kgyr_l200_bkg_mode(energy, mode, samples, partitions; kwargs...)
    m = mode
    @info "calculating posterior for energy $energy --- ", sum(partitions.exposure)
    f = s -> counts_density_obs(partitions, energy, s.S, m.B_l200a_all, s.αe_all; kwargs...)
    return broadcast(f, samples.v) ./ sum(partitions.exposure)
end

# this function calculates the counts density/(kg yr) at the best fit value
function counts_density_kgyr_l200_mode(energy, best_fit, partitions; kwargs...)
    n_cts = counts_density_obs(partitions, energy, best_fit.S, best_fit.B_l200a_all, best_fit.αe_all; kwargs...)
    return n_cts ./ sum(partitions.exposure)
end


function plot_l200_result(samples, part_event_index, events, partitions, fit_ranges, config)
    
    # we don't need all samples
    _samples = samples[1:20_000]
    _weights = weights(_samples.weight)

    # for the background, we just need one (random) energy point to compute the 68% smallest CI
    posterior = counts_density_kgyr_l200_post(2000, _samples, partitions, signal=false)
    _int = BAT.smallest_credible_intervals(posterior, _weights, nsigma_equivalent=1)
    length(_int) != 1 && @warn "[$E keV] 68% interval is disjoint"
    b_68 = first(_int)

    # let's use a narrow region around Qbb for the 90% CI on the signal
    energies = 2034:1.2:2043
    s_90 = Vector(undef, length(energies))

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
    b_mode = counts_density_kgyr_l200_mode(2000, best_fit_pars, partitions, signal=false)
    
    # let's speed this up with multi-threading
    Threads.@threads for i in 1:length(energies)
        E = energies[i]
        posterior = counts_density_kgyr_l200_bkg_mode(E, best_fit_pars, _samples, partitions)#, bkg=false)
        # no need for smallest interval here
        s_90[i] = b_mode..quantile(posterior, _weights, 0.9)
    end

    fig, ax = plt.subplots(figsize=(5, 2.5))

    # ...background 68%
    band = plt.fill_between(
        [1930, 2190], fill(b_68.left, 2), fill(b_68.right, 2),
        linewidth=0, alpha=0.5, color="#228833",
        label="68% C.I.",
        zorder=0,
    )
    # ...background best fit
    line, = plt.plot(
        [1930, 2190], [b_mode, b_mode],
        color="#228833",
        zorder=1,
    )

    ax.set_xlim(1930, 2190)
    ax.set_xlabel("Energy [keV]", loc="right")
    ax.set_ylabel("Counts / (keV kg yr)")

    # ...signal 90%
    limit_area = plt.fill_between(
        energies,
        getproperty.(s_90, :left), getproperty.(s_90, :right),
        linewidth=0, color="#4477AA",
        label="90% C.I.",
    )

    # plot the data
    unbinned = true

    # get energies of surviving events
    unflat_x = vcat([y for y in events if !isempty(y)])
    x = vcat(unflat_x...)
    data_art = nothing
    exp_tot = sum(partitions.exposure)

    if unbinned
        data_art, stemlines, baseline = ax.stem(
            x, fill(1e-3, length(x)), "#CC3311",
            basefmt="none",
        )
        data_art.set_markersize(1.5)
        # NOTE: comment the following if the stem has a length!
        data_art.set_zorder(99)
        data_art.set_clip_on(false)

        stemlines.set_linewidth(1)
        ax.set_ylim(1e-4, 2)
        ax.set_yscale("log")
    else
        w = 1/exp_tot
        _, _, data_art = plt.hist(
            x,
            weights=fill(w, length(x)),
            bins=1930:1:2190,
            color="#CC3311",
        )
        ax.set_ylim(1E-4, 2)
        ax.set_yscale("log")
    end

    # now this mess to plot the excluded regions
    plt.hist(
        [2104, 2119],
        bins=[2099, 2109, 2114, 2124],
        weights=[9, 9],
        color="black", alpha=0.4,
    )

    ax.legend(
        [data_art, (band, line), limit_area],
        [
            string(raw"LEGEND-200 data [", round(exp_tot, digits=2), " kg yr]"),
            raw"Background best fit and 68% C.I. interval",
            raw"$T^{0\nu}_{1/2} > 0.4 \times 10^{26}$ yr [90% C.I.]",
        ],
        loc="upper left",
        fontsize="small",
    )

    plt.tight_layout()
    plt.savefig("l200-result.pdf", bbox_inches="tight")
end
