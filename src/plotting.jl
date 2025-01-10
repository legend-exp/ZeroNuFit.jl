using Random, LinearAlgebra, Statistics, Distributions, StatsBase
using Plots
using BAT, DensityInterface, IntervalSets
using TypedTables
using Cuba
using Base.Filesystem
using PDFmerger: append_pdf!
using ColorSchemes
using OrderedCollections

default(
    framestyle = :box,               # Grid line transparency
    background_color = :white,       # Background color of the plot,
    titlefontsize = 10,     # Global title font size
    guidefontsize = 10,     # Global axis label font size
    tickfontsize = 10,      # Global tick label font size
    legendfontsize = 8,     # Global legend font size
)

tol_colors = ColorSchemes.tol_muted
color_schemes = Dict(
    :blue => [tol_colors[1], tol_colors[3], tol_colors[2]],
    :default => BAT.default_colors,
    :red => [:red4, :red, :salmon],
    :green => [:darkgreen, :chartreuse3, :palegreen1],
    :purple => [tol_colors[8], tol_colors[9], tol_colors[7]],
    :muted => [:olivedrab, :goldenrod, :indianred1],
)

function fit_model(
    idx_part_with_events::Int,
    part_k::NamedTuple,
    p::NamedTuple,
    settings::Dict,
    bkg_shape::Symbol,
    fit_range,
    x,
)
    Qbb = constants.Qbb
    N_A = constants.N_A
    m_76 = constants.m_76
    sig_units = constants.sig_units

    deltaE = sum([arr[2] - arr[1] for arr in fit_range])
    eff = nothing

    # logic for efficiency it can be either correlated, uncorrelated or fixed
    if settings[:bkg_only] == true
        eff = 0

    elseif (settings[:eff_correlated] == true)
        eff_group = part_k.eff_name
        eff = part_k.eff_tot + p[eff_group] * part_k.eff_tot_sigma

    elseif (
        idx_part_with_events != 0 &&
        settings[:eff_correlated] == false &&
        settings[:eff_fixed] == false
    )

        eff = p.ε[idx_part_with_events]
    else
        eff = part_k.eff_tot
    end

    b_name = part_k.bkg_name
    model_b_k = deltaE * part_k.exposure * p[b_name]
    term1 = model_b_k * get_bkg_pdf(bkg_shape, x, p, b_name, fit_range)

    if (settings[:bkg_only] == false)
        reso, bias = get_energy_scale_pars(part_k, p, settings, idx_part_with_events)
        model_s_k = log(2) * N_A * part_k.exposure * (eff) * (p.S * sig_units) / m_76
        term2 =
            model_s_k *
            get_signal_pdf(part_k.signal_name, x, Qbb, bias, reso, part_k, fit_range)
    else
        term2 = 0
    end

    return term1 + term2

end


"""
    plot_data(hist::Histogram,name,partitions,part_event_index,pars,samples,posterior,plotflag,settings::Dict,bkg_shape::Symbol,fit_ranges)

Function to plot events in the Qbb analysis window and BAT fit results
"""
function plot_data(
    hist::Histogram,
    name,
    partitions,
    part_event_index,
    pars,
    samples,
    posterior,
    plotflag,
    settings::Dict,
    bkg_shape::Symbol,
    fit_ranges,
)

    counts = sum(hist.weights)
    p = plot()

    ymax = 1.5 * maximum(hist.weights)
    bin_edges = hist.edges[1]

    # flatten all nested arrays into a single array of numbers
    all_values = reduce(vcat, reduce(vcat, values(fit_ranges)))
    # find the minimum and maximum values
    min_x = minimum(all_values)
    max_x = maximum(all_values)

    plot!(
        p,
        hist,
        st = :steps,
        label = "Data",
        title = "$name",
        xlabel = "Energy (keV)",
        ylabel = "Counts/ (1 keV)",
        ylim = (0, ymax),
        xlim = (min_x, max_x),
        color = "dark blue",
        fill = true,
        framestyle = :box,
    )

    #plot fit model
    function build_model_for_plotting(params, x)
        model = 0
        for (idx_k, part_k) in enumerate(partitions)
            fit_range = fit_ranges[part_k.fit_group]
            if part_event_index[idx_k] != 0
                idx_k_with_events = part_event_index[idx_k]
                model +=
                    fit_model(
                        idx_k_with_events,
                        part_k,
                        params,
                        settings,
                        bkg_shape,
                        fit_range,
                        x,
                    ) * diff(bin_edges)[1]
            else
                model +=
                    fit_model(0, part_k, params, settings, bkg_shape, fit_range, x) *
                    diff(bin_edges)[1]
            end
        end
        return model
    end

    if plotflag["bandfit_and_data"]
        plot!(
            p,
            min_x:0.1:max_x,
            build_model_for_plotting,
            samples,
            alpha = 0.4,
            median = false,
            globalmode = false,
            fillalpha = 0.3,
        ) #TO DO: take only some samples
        _, best_fit_pars = get_global_mode(samples, posterior)
        plot!(
            p,
            min_x:0.1:max_x,
            x -> build_model_for_plotting(best_fit_pars, x),
            label = "Fit",
            lw = 2,
            color = "red",
        )
    else
        _, best_fit_pars = get_global_mode(samples, posterior)
        plot!(
            p,
            min_x:0.1:max_x,
            x -> build_model_for_plotting(best_fit_pars, x),
            label = "Fit",
            lw = 2,
            color = "red",
        )
    end

    # exclude the gamma lines
    shape_x = [
        constants.gamma_2113_keV,
        constants.gamma_2113_keV,
        constants.gamma_2123_keV,
        constants.gamma_2123_keV,
    ]
    shape_x2 = [
        constants.gamma_2098_keV,
        constants.gamma_2098_keV,
        constants.gamma_2108_keV,
        constants.gamma_2108_keV,
    ]
    shape_x3 = [
        constants.gamma_2199_keV,
        constants.gamma_2199_keV,
        constants.gamma_2209_keV,
        constants.gamma_2209_keV,
    ]
    shape_y = [0, ymax, ymax, 0]

    plot!(
        p,
        shape_x,
        shape_y,
        fillalpha = 0.4,
        fill = true,
        line = false,
        linecolor = :transparent,
        fillcolor = :gray,
        label = "Excluded band",
    )
    plot!(
        p,
        shape_x2,
        shape_y,
        fillalpha = 0.4,
        fill = true,
        line = false,
        linecolor = :transparent,
        fillcolor = :gray,
        label = false,
    )
    plot!(
        p,
        shape_x3,
        shape_y,
        fillalpha = 0.4,
        fill = true,
        line = false,
        linecolor = :transparent,
        fillcolor = :gray,
        label = false,
    )

    return p
end



function plot_fit_and_data(
    partitions,
    events,
    part_event_index,
    samples,
    posterior,
    pars,
    output,
    config,
    fit_ranges;
    toy_idx = nothing,
)

    plotflag = config["plot"]
    settings = get_settings(config)
    bkg_shape, _ = get_bkg_info(config)

    # create histo with energies 
    energies = []
    for (idx_k, part_k) in enumerate(partitions)
        if events[idx_k] != Any[]
            for energy in events[idx_k]
                append!(energies, energy)
            end
        end
    end

    # find the minimum and maximum x values
    all_values = reduce(vcat, reduce(vcat, values(fit_ranges)))
    min_x = minimum(all_values)
    max_x = maximum(all_values)

    hist_data = append!(Histogram(min_x:1:max_x), energies)
    p_fit = plot_data(
        hist_data,
        "",
        partitions,
        part_event_index,
        pars,
        samples,
        posterior,
        plotflag,
        settings,
        bkg_shape,
        fit_ranges,
    )

    log_suffix = toy_idx == nothing ? "" : "_$(toy_idx)"
    savefig(joinpath(output, "plots/fit_over_data$log_suffix.pdf"))

end


"""
    plot_correlation_matrix(samples,output;par_names=nothing,toy_idx=nothing)

Plots the correlation matrixs
"""
function plot_correlation_matrix(samples, output; par_names = nothing, toy_idx = nothing)
    unshaped_samples, f_flatten = bat_transform(Vector, samples)
    covariance_matrix = cov(unshaped_samples)
    var = std(unshaped_samples)

    corr = 100 * sqrt(covariance_matrix ./ (var .* var'))
    heatmap(
        corr,
        xlabel = "Parameter Index",
        ylabel = "Parameter Index",
        color = :diverging_bwr_40_95_c42_n256,
        clim = (-100, 100),
        ctitle = "Correlation Coefficient",
    )

    log_suffix = toy_idx == nothing ? "" : "_$(toy_idx)"
    savefig(joinpath(output, "plots/correlations$log_suffix.pdf"))
end

function plot_two_dim_posteriors(
    samples,
    pars,
    output;
    par_names = nothing,
    toy_idx = nothing,
)
    first_sample = samples.v[1]

    for par_x in pars
        par_entry = first_sample[par_x]

        if (length(par_entry) != 1 || !(par_entry isa AbstractFloat))
            continue
        end
        for par_y in pars
            if (par_x == par_y)
                continue
            end
            par_entry = first_sample[par_y]
            if (length(par_entry) != 1 || !(par_entry isa AbstractFloat))
                continue
            end

            x = get_par_posterior(samples, par_x, idx = nothing)
            y = get_par_posterior(samples, par_y, idx = nothing)

            p = histogram2d(
                x,
                y,
                bins = 200,
                cmap = :batlow,
                xlabel = par_names[par_x],
                ylabel = par_names[par_y],
                right_margin = 10Plots.mm,
            )
            log_suffix = toy_idx == nothing ? "" : "_$(toy_idx)"
            savefig(p, "temp$log_suffix.pdf")
            append_pdf!(
                joinpath(output, "plots/2D_posterior$log_suffix.pdf"),
                "temp$log_suffix.pdf",
                cleanup = true,
            )
        end
    end


end


##############################################
##############################################
##############################################

"""
    plot_marginal_distr(partitions,samples,pars,output;sqrt_prior=false,priors=nothing,par_names=nothing,plot_config=nothing,s_max=nothing,hier=false,toy_idx=nothing) 

Function to plot 1D and 2D marginalized distributions (and priors)
"""
function plot_marginal_distr(
    partitions,
    samples,
    pars,
    output;
    sqrt_prior = false,
    priors = nothing,
    par_names = nothing,
    plot_config = nothing,
    s_max = nothing,
    hier = false,
    toy_idx = nothing,
)

    name = split(output, "output/")[end]
    first_sample = samples.v[1]
    unshaped_samples, f_flatten = bat_transform(Vector, samples)
    @debug "Unshaped samples:", bat_report(unshaped_samples)

    # remove old file
    log_suffix = toy_idx == nothing ? "" : "_$(toy_idx)"
    if isfile(joinpath(output, "plots/marg_posterior$log_suffix.pdf"))
        Filesystem.rm(joinpath(output, "plots/marg_posterior$log_suffix.pdf"), force = true)
    end


    # get a color scheme
    if plot_config != nothing && haskey(plot_config, "scheme")
        color_scheme = color_schemes[Symbol(plot_config["scheme"])]
    else
        color_scheme = BAT.default_colors
    end
    if plot_config != nothing && haskey(plot_config, "alpha")
        alpha = plot_config["alpha"]
    else
        alpha = 1
    end

    # 1D posteriors
    ct = 1
    for par in pars
        par_entry = first_sample[par]

        # checking if it is a 'AbstractFloat' helps avoiding cases were multivariate parameters have 1 entry only
        if (length(par_entry) == 1 && par_entry isa AbstractFloat)

            post = get_par_posterior(samples, par, idx = nothing)
            if (par == :S || par == :B)
                mini = 0
            else
                mini = minimum(post)
            end

            if (par_names != nothing)
                xname = par_names[par]
            end

            p = plot(
                samples,
                par,
                mean = false,
                std = false,
                globalmode = true,
                marginalmode = true,
                nbins = 200,
                xlim = (mini, maximum(post)),
                colors = color_scheme,
                alpha = alpha,
                lw = 0,
                linecolor = :black,
            )
            xaxis!(xname)
            yaxis!("Probability Density")
            ylims!(0, ylims()[2])

            x = range(mini, stop = maximum(post), length = 1000)
            maxi = ylims()[2]
            # plot prior
            if priors != nothing
                if (par == :S && sqrt_prior)
                    y = x -> 1 ./ (2 * sqrt(s_max * x))
                else
                    color = "black"
                    if (hier == false)
                        y = pdf(priors[par], x)
                        plot!(x, y, label = "prior", color = "black")

                    elseif (haskey(priors.pdist, par))

                        y = pdf(priors.pdist[par], x)
                        plot!(x, y, label = "prior", color = "black")

                    else
                        for i = 1:50
                            rando = rand(priors.pdist)
                            distribution = Categorical(samples.weight / sum(samples.weight))
                            random_index = rand(distribution)
                            rando = samples.v[random_index]
                            y = pdf(priors.f(rando)[par], x)
                            maxi = maximum([maximum(y), maxi])
                            color = "grey"

                            plot!(x, y, color = "grey", alpha = 0.3, label = nothing)
                            ylims!(0, maxi)

                        end
                    end
                end

            end

            log_suffix = toy_idx == nothing ? "" : "_$(toy_idx)"
            savefig(p, "temp$log_suffix.pdf")
            append_pdf!(
                joinpath(output, "plots/marg_posterior$log_suffix.pdf"),
                "temp$log_suffix.pdf",
                cleanup = true,
            )
            ct += 1

            # multivariate parameters    
        else
            for idx = 1:length(par_entry)
                post = get_par_posterior(samples, par, idx = idx)

                xlab = string("$(par)[$(idx)]")
                ylab = string("Probability Density")
                if (par_names != nothing)
                    xname = par_names[par][idx]
                end

                p = plot(
                    unshaped_samples,
                    ct,
                    mean = false,
                    std = false,
                    globalmode = true,
                    marginalmode = true,
                    nbins = 200,
                    xlabel = xlab,
                    ylabel = ylab,
                    xlim = (minimum(post), maximum(post)),
                    colors = color_scheme,
                    alpha = alpha,
                    linecolor = :black,
                )
                mini = minimum(post)
                x = range(mini, stop = maximum(post), length = 1000)

                if priors != nothing

                    if (hier == true && haskey(priors.pdist, par))
                        y = pdf(priors.pdist[par].v[idx], x)
                    elseif (hier == true)
                        rando = rand(priors.pdist)
                        rando = samples.v[1]
                        y = pdf(priors.f(rando)[par].v[idx], x)
                    else
                        y = pdf(priors[par].v[idx], x)
                    end
                    plot!(x, y, label = "prior", color = "black")
                end

                xaxis!(xname)
                ylims!(0, ylims()[2])
                log_suffix = toy_idx == nothing ? "" : "_$(toy_idx)"
                savefig(p, "temp$log_suffix.pdf")
                append_pdf!(
                    joinpath(output, "plots/marg_posterior$log_suffix.pdf"),
                    "temp$log_suffix.pdf",
                    cleanup = true,
                )
                ct += 1
            end
        end
    end

end
