### mbb.jl
#
# Authors: Sofia Calgaro, Toby Dixon
# 
###
using Pkg
Pkg.activate(".")
Pkg.instantiate()
using ArgParse
using JSON
using DataStructures
using PropDicts
using BAT
import HDF5
using TypedTables
using Plots, LaTeXStrings
using Base.Filesystem
using PDFmerger: append_pdf!
using Plots
using ColorSchemes
using Random, LinearAlgebra, Statistics, Distributions, StatsBase

# style
default(
    framestyle = :box,
    background_color = :white,
    titlefontsize = 10,
    guidefontsize = 10,
    tickfontsize = 10,
    legendfontsize = 10,
    linewidth = 1.5,
)
tol_colors = ColorSchemes.tol_muted
tol_mk15 = ColorSchemes.mk_15

const nme_belley_central = 2.60 # A. Belley et al., Phys. Rev. Lett. 132, 182502 (2024), arXiv:2308.15634
const nme_belley_up = 1.28      # A. Belley et al., Phys. Rev. Lett. 132, 182502 (2024), arXiv:2308.15634
const nme_belley_low = 1.36     # A. Belley et al., Phys. Rev. Lett. 132, 182502 (2024), arXiv:2308.15634
const nme_adams_up = 6.34       # C. Adams et al., (2022), arXiv:2212.11099
const nme_adams_low = 2.35      # C. Adams et al., (2022), arXiv:2212.11099
const nme_gerda_up = 6.04       # GERDA Phys. Rev. Lett. 125 (2020) 252502 [arXiv:2009.06079]
const nme_gerda_low = 2.66      # GERDA Phys. Rev. Lett. 125 (2020) 252502 [arXiv:2009.06079]
const phase_space = 0.237 * 10^-14
const gA = 1.2724

# inputs
nme_belley = (central = nme_belley_central, up = nme_belley_up, low = nme_belley_low)
nme_adams = (up = nme_adams_up, low = nme_adams_low)
nme_gerda = (up = nme_gerda_up, low = nme_gerda_low)
phase_space = phase_space
gA = gA

function mbb(S; G = 1, M = 1, gA = gA, me = me_keV * 1000)
    return sqrt.(S) .* M .^ -1 * sqrt(10.0^-27 * me^2 * G^-1 * gA^-4)
end


function mbb2(S; G = 1, M = 1, gA = gA, me = me_keV * 1000)
    return S .* M .^ -2 * (10.0^-27 * me^2 * G^-1 * gA^-4)
end

function get_S_posterior(samples)
    S_samples = []

    for samp in samples
        v = samp.v
        weight = samp.weight
        for w = 1:1:weight
            append!(S_samples, v[:S])
        end
    end
    S_samples
end

function get_asymmetric_gaussian(S_samples, xmin, xmax, mu, sigma_low, sigma_high)
    sample = rand(Uniform(xmin, xmax), length(S_samples))
    asymm_gaussian_samples = Float64[]
    for el in sample
        if rand() < mu
            y = rand(Truncated(Normal(mu, sigma_low), xmin, xmax))
            push!(asymm_gaussian_samples, y)
        else
            y = rand(Truncated(Normal(mu, sigma_high), xmin, xmax))
            push!(asymm_gaussian_samples, y)
        end
    end
    return asymm_gaussian_samples
end

using DelimitedFiles, StatsBase, Distributions, Plots
function plot_true_Belley()

    # Load the data, skipping the first 3 rows and using the second column (index 2)
    data = readdlm("attic/belley_posterior_PDF_samples.csv", ',', skipstart = 3)[:, 2]
    data = data[data.>=0]
    bins = range(0, 10, length = 200)
    binwidth = step(bins)

    plot_posterior(
        string("M") * L"_{0\nu}",
        data,
        10,
        string("M") * L"_{0\nu}",
        "",
        false,
        ".",
        "true_belley.pdf",
        line = nme_belley.central,
    )
    plot_posterior(
        string("M") * L"_{0\nu}^{-1}",
        data .^ -1,
        1,
        string("M") * L"_{0\nu}^{-1}",
        "",
        false,
        ".",
        "true_belley.pdf",
        line = 1 / nme_belley.central,
    )
    plot_posterior(
        string("M") * L"_{0\nu}^{-2}",
        data .^ -2,
        1,
        string("M") * L"_{0\nu}^{-2}",
        "",
        false,
        ".",
        "true_belley.pdf",
        line = 1 / nme_belley.central / nme_belley.central,
    )

    # sampling from Belley's NME posterior
    M_dist =
        Truncated(Normal(nme_belley.central, (nme_belley.up + nme_belley.low) / 2), 0, 30)
    M_samples = rand(M_dist, length(data))
    p = plot()
    vline!([2.60], linestyle = :dash, color = :gray, label = "Belley's central value")
    hist = append!(Histogram(0:10/200:10), data)
    plot!(
        hist,
        st = :steps,
        label = "Belley's pdf",
        color = tol_colors[1],
        linewidth = 1.8,
        tight = true,
    )
    hist = append!(Histogram(0:10/200:10), M_samples)
    plot!(
        hist,
        st = :steps,
        label = "Symmetric Gaussian",
        color = tol_colors[2],
        linewidth = 1.5,
        tight = true,
    )

    xaxis!(String("M") * L"_{0\nu}")
    yaxis!("Proability Density")
    xlims!(0, 10)
    ylims!(0, ylims()[2])
    savefig(p, "belley_vs_symm_gaussian.pdf")

end

function plot_posterior(
    name::Union{String,LaTeXString},
    post,
    max_x,
    axis_label::Union{String,LaTeXString},
    label::String,
    quantiles::Bool,
    output_path::String,
    filename::String;
    line = nothing,
    post_other = nothing,
    label_other = nothing,
    col1 = tol_colors[1],
    col2 = tol_colors[2],
    col3 = tol_colors[3],
    col4 = tol_colors[4],
)

    # plot the first histogram
    hist = append!(Histogram(0:max_x/200:max_x), post)
    p = plot(
        hist,
        st = :steps,
        title = "Posterior distribution of $name",
        label = label,
        color = col1,
        linewidth = 1.5,
        tight = true,
    )
    if (quantiles)
        vline!(
            [quantile(post, 0.9)],
            linestyle = :dash,
            color = col1,
            label = "90% quantile",
        )
    end
    if (line != nothing)
        vline!([line], linestyle = :dash, color = col4, label = "Central NME value")
    end
    # plot the second histogram
    if (post_other != nothing)
        hist_other = append!(Histogram(0:max_x/200:max_x), post_other)

        plot!(
            hist_other,
            st = :steps,
            title = "Posterior distribution",
            label = label_other,
            color = "orange",
        )
        if (quantiles)
            vline!(
                [quantile(post_other, 0.9)],
                linestyle = :dash,
                color = "orange",
                label = "90% quantile",
            )
        end
    end
    xaxis!(axis_label)
    yaxis!("Proability Density")
    xlims!(0, max_x)
    ylims!(0, ylims()[2])

    savefig(p, "temporary_pdf.pdf")
    append_pdf!(joinpath(output, filename), "temporary_pdf.pdf", cleanup = true)

end


function plot_mbb_plots(file_path::String, output_path::String; filename = "mbb.pdf")

    # remove previous mbb pdf if already present
    if isfile(joinpath(output_path, filename))
        Filesystem.rm(joinpath(output_path, filename), force = true)
    end

    # get the relevant samples
    samples = bat_read(file_path).result
    S_samples = get_S_posterior(samples)

    # sampling from Belley's NME posterior
    data = readdlm("inputs/belley_posterior_PDF_samples.csv", ',', skipstart = 3)[:, 2]
    data = data[data.>=0]
    M_samples = data
    # outdated symm truncated gaussian
    #M_dist = Truncated(Normal(nme_belley.central, (nme_belley.up + nme_belley.low) / 2), 0, 30)
    #M_samples = rand(M_dist, length(S_samples))

    # asymmetric Belley's pdf
    M_asymm_samples = get_asymmetric_gaussian(
        S_samples,
        0,
        30,
        nme_belley.central,
        nme_belley.low,
        nme_belley.up,
    )

    # uniform sampling from min-max NME values (GERDA ones)
    NME_dist = Uniform(nme_gerda.low, nme_gerda.up)
    NME_samples = rand(NME_dist, length(S_samples))

    # fixed NME (Belley's values)
    fixed_belley_mbb_samples = mbb(S_samples, G = phase_space, M = nme_belley.central)
    fixed_belley_mbb_samples_low =
        mbb(S_samples, G = phase_space, M = nme_belley.central - nme_belley.low)
    fixed_belley_mbb_samples_up =
        mbb(S_samples, G = phase_space, M = nme_belley.central + nme_belley.up)
    fixed_belley_mbb2_samples = mbb2(S_samples, G = phase_space, M = nme_belley.central)
    # fixed NME (Adam's values)
    fixed_adams_mbb_samples_up = mbb(S_samples, G = phase_space, M = nme_adams.up)
    fixed_adams_mbb_samples_low = mbb(S_samples, G = phase_space, M = nme_adams.low)
    fixed_adams_mbb2_samples_up = mbb2(S_samples, G = phase_space, M = nme_adams.up)
    fixed_adams_mbb2_samples_low = mbb2(S_samples, G = phase_space, M = nme_adams.low)
    # fixed NME (GERDA final results)
    fixed_gerda_mbb_samples_up = mbb(S_samples, G = phase_space, M = nme_gerda.up)
    fixed_gerda_mbb_samples_low = mbb(S_samples, G = phase_space, M = nme_gerda.low)
    fixed_gerda_mbb2_samples_up = mbb2(S_samples, G = phase_space, M = nme_gerda.up)
    fixed_gerda_mbb2_samples_low = mbb2(S_samples, G = phase_space, M = nme_gerda.low)

    # vary the NME according to Belley's posterior
    marginalised_mbb_samples_Belley = mbb(S_samples, G = phase_space, M = M_samples)
    marginalised_mbb2_samples_Belley = mbb2(S_samples, G = phase_space, M = M_samples)
    marginalised_mbb_samples_Belley_asymm =
        mbb(S_samples, G = phase_space, M = M_asymm_samples)
    marginalised_mbb2_samples_Belley_asymm =
        mbb2(S_samples, G = phase_space, M = M_asymm_samples)
    # vary the NME according to uniform sampled NMes
    marginalised_mbb_samples_uniformNME = mbb(S_samples, G = phase_space, M = NME_samples)
    marginalised_mbb2_samples_uniformNME = mbb2(S_samples, G = phase_space, M = NME_samples)

    @info "90% mbb quantile, fixing NME to Belley's best value: ",
    quantile(fixed_belley_mbb_samples, 0.9) * 1000
    @info "90% mbb quantile, fixing NME to Belley's best value + 1 sigma: ",
    quantile(fixed_belley_mbb_samples_up, 0.9) * 1000
    @info "90% mbb quantile, fixing NME to Belley's best value - 1 sigma: ",
    quantile(fixed_belley_mbb_samples_low, 0.9) * 1000
    @info "90% mbb quantile, varying Belley's NME (symm): ",
    quantile(marginalised_mbb_samples_Belley, 0.9) * 1000
    @info "90% mbb quantile, varying Belley's NME (asymm): ",
    quantile(marginalised_mbb_samples_Belley_asymm, 0.9) * 1000
    @info "90% mbb quantile, fixing NME to lowest GERDA's NME value: ",
    quantile(fixed_gerda_mbb_samples_low, 0.9) * 1000
    @info "90% mbb quantile, fixing NME to highest GERDA's NME value: ",
    quantile(fixed_gerda_mbb_samples_up, 0.9) * 1000
    @info "90% mbb quantile, fixing NME to lowest ADAM's NME value: ",
    quantile(fixed_adam_mbb_samples_low, 0.9) * 1000
    @info "90% mbb quantile, fixing NME to highest ADAM's NME value: ",
    quantile(fixed_adam_mbb_samples_up, 0.9) * 1000
    @info "90% mbb quantile, varying NME uniformly between min and max values: ",
    quantile(marginalised_mbb_samples_uniformNME, 0.9) * 1000

    # fix S and vary NME
    fix_S_samples = fill(20, 60000)
    M_samples = rand(M_dist, length(fix_S_samples))
    fixed_S_marginalised_mbb_samples = mbb(fix_S_samples, G = phase_space, M = M_samples)
    fixed_S_marginalised_mbb2_samples = mbb2(fix_S_samples, G = phase_space, M = M_samples)

    fixed_S_mbb_samples = mbb(fix_S_samples, G = phase_space, M = nme_belley.central)[1]
    fixed_S_mbb2_samples = mbb2(fix_S_samples, G = phase_space, M = nme_belley.central)[1]

    # SIGNAL posteriors
    plot_posterior(
        "S",
        S_samples,
        40,
        string("S [") * L"10^{-27} yr^{-1}]",
        "",
        true,
        output_path,
        filename,
    )
    plot_posterior(
        L"\sqrt{S}",
        S_samples .^ 0.5,
        10,
        L"\sqrt{S}" * string(" [") * L"10^{-27} yr^{-1/2}]",
        "",
        true,
        output_path,
        filename,
    )

    # NME posteriors
    plot_posterior(
        string("M") * L"_{0\nu} (Belley)",
        M_samples,
        10,
        string("M") * L"_{0\nu}",
        "",
        false,
        output_path,
        filename,
        line = nme_belley.central,
    )
    plot_posterior(
        string("M") * L"_{0\nu}^{-1} (Belley)",
        M_samples .^ -1,
        1,
        string("M") * L"_{0\nu}^{-1}",
        "",
        false,
        output_path,
        filename,
        line = 1 / nme_belley.central,
    )
    plot_posterior(
        string("M") * L"_{0\nu}^{-2} (Belley)",
        M_samples .^ -2,
        1,
        string("M") * L"_{0\nu}^{-2}",
        "",
        false,
        output_path,
        filename,
        line = 1 / (nme_belley.central^2),
    )
    # uniformly sampled NME
    plot_posterior(
        string("M") * L"_{0\nu} (uniform)",
        NME_samples,
        nme_adams.up + 1,
        string("M") * L"_{0\nu}",
        "",
        false,
        output_path,
        filename,
        col1 = tol_colors[2],
    )
    plot_posterior(
        string("M") * L"_{0\nu}^{-1} (uniform)",
        NME_samples .^ -1,
        (nme_adams.low - 1) .^ -1,
        string("M") * L"_{0\nu}^{-1}",
        "",
        false,
        output_path,
        filename,
        col1 = tol_colors[2],
    )
    plot_posterior(
        string("M") * L"_{0\nu}^{-2} (uniform)",
        NME_samples .^ -2,
        (nme_adams.low - 1) .^ -2,
        string("M") * L"_{0\nu}^{-2}",
        "",
        false,
        output_path,
        filename,
        col1 = tol_colors[2],
    )

    # Fixed NME - mbb
    plot_posterior(
        L"m_{\beta\beta}" * String(" "),
        fixed_belley_mbb_samples,
        1.0,
        string("m") * L"_{\beta\beta}" * string(" [eV]"),
        "Fixed NME",
        true,
        output_path,
        filename,
        post_other = marginalised_mbb_samples_Belley,
        label_other = "Marginalised NME",
    )

    # Fixed NME - mbb^2
    plot_posterior(
        L"m_{\beta\beta}^2" * String(" "),
        fixed_belley_mbb2_samples,
        0.4,
        string("m") * L"_{\beta\beta}^2" * string(" [eV") * L"^2]",
        "Fixed NME",
        true,
        output_path,
        filename,
        post_other = marginalised_mbb2_samples_Belley,
        label_other = "Marginalised NME",
    )

    # Fixed NME - mbb - uniform NME sampling VS max NME value
    plot_posterior(
        L"m_{\beta\beta}" * String(" "),
        marginalised_mbb_samples_uniformNME,
        0.7,
        string("m") * L"_{\beta\beta}" * string(" [eV]"),
        "Uniformly distributed NME",
        true,
        output_path,
        filename,
        col1 = tol_colors[2],
        post_other = fixed_adam_mbb_samples_up,
        label_other = "NME fixed to max. value " * string(nme_adams.up),
    )
    # Fixed NME - mbb - uniform NME sampling VS min NME value
    plot_posterior(
        L"m_{\beta\beta}" * String(" "),
        marginalised_mbb_samples_uniformNME,
        0.7,
        string("m") * L"_{\beta\beta}" * string(" [eV]"),
        "Uniformly distributed NME",
        true,
        output_path,
        filename,
        col1 = tol_colors[2],
        post_other = fixed_adam_mbb_samples_low,
        label_other = "NME fixed to min. value " * string(nme_adams.low),
    )

    # Fixed NME - mbb^2 - uniform NME sampling VS max NME value
    plot_posterior(
        L"m_{\beta\beta}^2" * String(" "),
        marginalised_mbb2_samples_uniformNME,
        0.1,
        string("m") * L"_{\beta\beta}^2" * string(" [eV") * L"^2]",
        "Uniformly distributed NME",
        true,
        output_path,
        filename,
        col1 = tol_colors[2],
        post_other = fixed_adam_mbb2_samples_up,
        label_other = "NME fixed to max. value " * string(nme_adams.up),
    )
    # Fixed NME - mbb^2 - uniform NME sampling VS min NME value
    plot_posterior(
        L"m_{\beta\beta}^2" * String(" "),
        marginalised_mbb2_samples_uniformNME,
        0.1,
        string("m") * L"_{\beta\beta}^2" * string(" [eV") * L"^2]",
        "Uniformly distributed NME",
        true,
        output_path,
        filename,
        col1 = tol_colors[2],
        post_other = fixed_adam_mbb2_samples_low,
        label_other = "NME fixed to min. value " * string(nme_adams.low),
    )

    # Fixed S - mbb
    plot_posterior(
        L"m_{\beta\beta}" * String(" (S=20 fixed)"),
        fixed_S_marginalised_mbb_samples,
        1.5,
        string("m") * L"_{\beta\beta}" * string(" [eV]"),
        "Marginalised NME",
        true,
        output_path,
        filename,
        line = fixed_S_mbb_samples,
        post_other = nothing,
        label_other = nothing,
    )

    # Fixed S - mbb^2
    plot_posterior(
        L"m_{\beta\beta}^2" * String(" (S=20 fixed)"),
        fixed_S_marginalised_mbb2_samples,
        0.7,
        string("m") * L"_{\beta\beta}^2" * string(" [eV") * L"^2]",
        "Marginalised NME",
        true,
        output_path,
        filename,
        line = fixed_S_mbb2_samples,
        post_other = nothing,
        label_other = nothing,
    )

end


function plot_mbb_Belley_studies(
    file_path::String,
    output_path::String;
    filename = "mbb.pdf",
)

    # remove previous mbb pdf if already present
    if isfile(joinpath(output_path, filename))
        Filesystem.rm(joinpath(output_path, filename), force = true)
    end

    # get the relevant samples
    samples = bat_read(file_path).result
    S_samples = get_S_posterior(samples)

    max_x = 0.8
    quantiles = false


    # m pdf for shifted uniformly sampled NMEs (fixed range width)
    p = plot()
    mycol = "orange"

    for (idx, nme) in enumerate([-2.5, -1.5, 0, 1.5, 2.5])
        maxNME = nme_adams.up + nme
        minNME = round(nme_adams.low + nme; digits = 2)
        NME_dist = Uniform(minNME, maxNME)
        uniform_samples = rand(NME_dist, length(S_samples))

        hist = append!(Histogram(0:10/200:10), uniform_samples)
        plot!(
            hist,
            st = :steps,
            title = String("M") * L"_{0\nu} (uniform)",
            label = "NME in $minNME-$maxNME",
            legend = :bottomright,
            linewidth = 1.5,
            fillalpha = 0.5,
            tight = true,
        )

    end

    xaxis!(String("M") * L"_{0\nu}")
    yaxis!("Proability Density")
    xlims!(0, 10)
    ylims!(0, ylims()[2])

    savefig(p, "temporary_pdf.pdf")
    append_pdf!(joinpath(output, filename), "temporary_pdf.pdf", cleanup = true)

    p = plot()
    mycol = "orange"

    for (idx, nme) in enumerate([-2.5, -1.5, 0, 1.5, 2.5])
        maxNME = nme_adams.up + nme
        minNME = round(nme_adams.low + nme; digits = 2)
        NME_dist = Uniform(minNME, maxNME)
        uniform_samples = rand(NME_dist, length(S_samples))

        hist = append!(Histogram(0:1/200:1), uniform_samples .^ -1)
        plot!(
            hist,
            st = :steps,
            title = String("M") * L"_{0\nu}^{-1} (uniform)",
            label = "NME in $minNME-$maxNME",
            legend = :topright,
            linewidth = 1.5,
            tight = true,
            fillalpha = 0.5,
        )

    end

    xaxis!(String("M") * L"_{0\nu}^{-1}")
    yaxis!("Proability Density")
    xlims!(0, 1)
    ylims!(0, ylims()[2])

    savefig(p, "temporary_pdf.pdf")
    append_pdf!(joinpath(output, filename), "temporary_pdf.pdf", cleanup = true)

    p = plot()
    for (idx, nme) in enumerate([-2.5, -1.5, 0, 1.5, 2.5])
        maxNME = nme_adams.up + nme
        minNME = round(nme_adams.low + nme; digits = 2)
        NME_dist = Uniform(minNME, maxNME)
        uniform_samples = rand(NME_dist, length(S_samples))

        hist = append!(Histogram(0:0.4/200:0.4), uniform_samples .^ -2)
        plot!(
            hist,
            st = :steps,
            title = String("M") * L"_{0\nu}^{-2} (uniform)",
            label = "NME in $minNME-$maxNME",
            legend = :topright,
            linewidth = 1.5,
            tight = true,
            fillalpha = 0.5,
        )

    end

    xaxis!(String("M") * L"_{0\nu}^{-2}")
    yaxis!("Proability Density")
    xlims!(0, 0.4)
    ylims!(0, ylims()[2])

    savefig(p, "temporary_pdf.pdf")
    append_pdf!(joinpath(output, filename), "temporary_pdf.pdf", cleanup = true)


    for (idx, nme) in enumerate([-2.5, -1.5, 0, 1.5, 2.5])
        p = plot()
        maxNME = nme_adams.up + nme
        minNME = round(nme_adams.low + nme; digits = 2)
        NME_dist = Uniform(minNME, maxNME)
        uniform_samples = rand(NME_dist, length(S_samples))

        hist = append!(Histogram(0:1/200:1), uniform_samples .^ -1)
        plot!(
            hist,
            st = :steps,
            title = String("M") * L"_{0\nu}^{-1} (uniform)",
            label = "NME in $minNME-$maxNME",
            legend = :topright,
            linewidth = 1.5,
            tight = true,
            fillalpha = 0.5,
        )
        xaxis!(String("M") * L"_{0\nu}^{-1}")
        yaxis!("Proability Density")
        xlims!(0, 1)
        ylims!(0, ylims()[2])

        savefig(p, "temporary_pdf.pdf")
        append_pdf!(joinpath(output, filename), "temporary_pdf.pdf", cleanup = true)
    end
    for (idx, nme) in enumerate([-2.5, -1.5, 0, 1.5, 2.5])
        p = plot()
        maxNME = nme_adams.up + nme
        minNME = round(nme_adams.low + nme; digits = 2)
        NME_dist = Uniform(minNME, maxNME)
        uniform_samples = rand(NME_dist, length(S_samples))

        hist = append!(Histogram(0:0.4/200:0.4), uniform_samples .^ -2)
        plot!(
            hist,
            st = :steps,
            title = String("M") * L"_{0\nu}^{-2} (uniform)",
            label = "NME in $minNME-$maxNME",
            legend = :topright,
            linewidth = 1.5,
            tight = true,
            fillalpha = 0.5,
        )
        xaxis!(String("M") * L"_{0\nu}^{-2}")
        yaxis!("Proability Density")
        xlims!(0, 0.4)
        ylims!(0, ylims()[2])

        savefig(p, "temporary_pdf.pdf")
        append_pdf!(joinpath(output, filename), "temporary_pdf.pdf", cleanup = true)
    end


    # m pdf for shifted uniformly sampled NMEs (fixed range width)
    p = plot()
    uniform_limits = []
    mycol = "orange"

    for (idx, nme) in enumerate([-2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5])
        maxNME = nme_adams.up + nme
        minNME = round(nme_adams.low + nme; digits = 2)
        NME_dist = Uniform(minNME, maxNME)
        uniform_samples = rand(NME_dist, length(S_samples))
        marginalised_mbb_samples_uniformNME =
            mbb(S_samples, G = phase_space, M = uniform_samples)
        append!(uniform_limits, quantile(marginalised_mbb_samples_uniformNME, 0.9))

        hist = append!(Histogram(0:0.3/200:0.3), marginalised_mbb_samples_uniformNME)
        plot!(
            hist,
            st = :steps,
            title = String("Posterior distribution of ") *
                    L"m_{\beta\beta}" *
                    String(" varying NME (uniform sampling)"),
            label = "NME in $minNME-$maxNME",
            #color = tol_mk15[idx],
            linewidth = 1.5,
            tight = true,
        )
        if (quantiles)
            vline!(
                [quantile(marginalised_mbb_samples_uniformNME, 0.9)],
                linestyle = :dash,
                #color = tol_mk15[idx],
                #label = "90% quantile",
            )
        end

    end

    xaxis!(L"m_{\beta\beta}" * String(" [eV]"))
    yaxis!("Proability Density")
    xlims!(0, 0.3)
    ylims!(0, ylims()[2])

    savefig(p, "temporary_pdf.pdf")
    append_pdf!(joinpath(output, filename), "temporary_pdf.pdf", cleanup = true)


    # m^2 pdf for shifted uniformly sampled NMEs (fixed range width)
    p = plot()
    mycol = "orange"

    for (idx, nme) in enumerate([-2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5])
        maxNME = nme_adams.up + nme
        minNME = round(nme_adams.low + nme; digits = 2)
        NME_dist = Uniform(minNME, maxNME)
        uniform_samples = rand(NME_dist, length(S_samples))
        marginalised_mbb2_samples_uniformNME =
            mbb2(S_samples, G = phase_space, M = uniform_samples)

        hist = append!(Histogram(0:0.02/200:0.02), marginalised_mbb2_samples_uniformNME)
        plot!(
            hist,
            st = :steps,
            title = String("Posterior distribution of ") *
                    L"m_{\beta\beta}^2" *
                    String(" varying NME (uniform sampling)"),
            label = "NME in $minNME-$maxNME",
            linewidth = 1.5,
            tight = true,
        )
        if (quantiles)
            vline!([quantile(marginalised_mbb2_samples_uniformNME, 0.9)], linestyle = :dash)
        end

    end

    xaxis!(L"m_{\beta\beta}^2" * String(" [eV]"))
    yaxis!("Proability Density")
    xlims!(0, 0.02)
    ylims!(0, ylims()[2])

    savefig(p, "temporary_pdf.pdf")
    append_pdf!(joinpath(output, filename), "temporary_pdf.pdf", cleanup = true)


    # limits vs NME for shifted uniformly sampled NMEs (fixed range width)
    p = plot()
    vline!([nme_adams.low], linestyle = :dash, color = "black")
    vline!([nme_adams.up], linestyle = :dash, color = "black")

    gerda_T12 = 1.83
    gerda_low_limit = mbb(1 / gerda_T12 * 10, G = phase_space, M = nme_adams_up) * 1000
    gerda_upp_limit = mbb(1 / gerda_T12 * 10, G = phase_space, M = nme_adams_low) * 1000
    hspan!(
        [gerda_low_limit, gerda_upp_limit],
        fillrange = 1,
        fillalpha = 0.2,
        color = :green,
    )

    for (idx, nme) in enumerate([-2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5])
        maxNME = nme_adams.up + nme
        minNME = round(nme_adams.low + nme; digits = 2)
        mycol = "orange"
        mycol2 = "navy"
        limit = uniform_limits[idx] * 1000
        plot!([minNME, maxNME], [limit, limit], linewidth = 2, legend = :false)
    end

    yaxis!(L"m_{\beta\beta}" * String(" upper limit, 90% CI [meV]"))
    xaxis!("NME value")

    savefig(p, "temporary_pdf.pdf")
    append_pdf!(joinpath(output, filename), "temporary_pdf.pdf", cleanup = true)



    # m pdf for gaussian sampled NMEs
    p = plot()
    limits = []
    mycol = "orange"

    for (idx, nme) in
        enumerate([0.5, 1.0, 1.3, 1.7, 2.0, nme_belley.central, 3.0, 4.0, 5.0, 7.5])
        M_dist = Truncated(Normal(nme, (nme_belley.up + nme_belley.low) / 2), 0, 30)
        M_samples = rand(M_dist, length(S_samples))
        marginalised_mbb_samples_Belley = mbb(S_samples, G = phase_space, M = M_samples)
        append!(limits, quantile(marginalised_mbb_samples_Belley, 0.9))

        hist = append!(Histogram(0:max_x/200:max_x), marginalised_mbb_samples_Belley)
        plot!(
            hist,
            st = :steps,
            title = String("Posterior distribution of ") *
                    L"m_{\beta\beta}" *
                    String(" varying NME (gaussian sampling)"),
            label = "NME = $nme",
            #color = tol_mk15[idx],
            linewidth = 1.5,
            tight = true,
        )
        if (quantiles)
            vline!(
                [quantile(marginalised_mbb_samples_Belley, 0.9)],
                linestyle = :dash,
                #color = tol_mk15[idx],
                #label = "90% quantile",
            )
        end

    end

    xaxis!(L"m_{\beta\beta}" * String(" [eV]"))
    yaxis!("Proability Density")
    xlims!(0, max_x)
    ylims!(0, ylims()[2])

    savefig(p, "temporary_pdf.pdf")
    append_pdf!(joinpath(output, filename), "temporary_pdf.pdf", cleanup = true)


    # m pdf for fixed NMEs
    p = plot()
    fixed_limits = []
    for (idx, nme) in
        enumerate([0.5, 1.0, 1.3, 1.7, 2.0, nme_belley.central, 3.0, 4.0, 5.0, 7.5])
        marginalised_mbb_samples_Belley = mbb(S_samples, G = phase_space, M = nme)
        append!(fixed_limits, quantile(marginalised_mbb_samples_Belley, 0.9))

        hist = append!(Histogram(0:max_x/200:max_x), marginalised_mbb_samples_Belley)
        plot!(
            hist,
            st = :steps,
            title = String("Posterior distribution of ") *
                    L"m_{\beta\beta}" *
                    String(" varying NME (fixed)"),
            label = "NME = $nme",
            #color = tol_mk15[idx],
            linewidth = 1.5,
            tight = true,
        )
        if (quantiles)
            vline!(
                [quantile(marginalised_mbb_samples_Belley, 0.9)],
                linestyle = :dash,
                #color = tol_mk15[idx],
                #label = "90% quantile",
            )
        end

    end

    xaxis!(L"m_{\beta\beta}" * String(" [eV]"))
    yaxis!("Proability Density")
    xlims!(0, max_x)
    ylims!(0, ylims()[2])

    savefig(p, "temporary_pdf.pdf")
    append_pdf!(joinpath(output, filename), "temporary_pdf.pdf", cleanup = true)

    fixed_68_low_limits = []
    fixed_68_upp_limits = []
    for (idx, nme) in
        enumerate([0.5, 1.0, 1.3, 1.7, 2.0, nme_belley.central, 3.0, 4.0, 5.0, 7.5])
        nme_low = nme - (nme_belley.up + nme_belley.low) / 2
        marginalised_mbb_low = mbb(S_samples, G = phase_space, M = nme_low)
        nme_up = nme + (nme_belley.up + nme_belley.low) / 2
        marginalised_mbb_up = mbb(S_samples, G = phase_space, M = nme_up)
        append!(fixed_68_low_limits, quantile(marginalised_mbb_low, 0.9))
        append!(fixed_68_upp_limits, quantile(marginalised_mbb_up, 0.9))
    end



    # M pdfs for gaussian sampled NMEs
    for (idx, nme) in enumerate([1.0, 1.5, 2.0, nme_belley.central, 4.0, 7.5])
        M_dist = Truncated(Normal(nme, (nme_belley.up + nme_belley.low) / 2), 0, 30)
        M_samples = rand(M_dist, length(S_samples))

        p = plot()
        hist = append!(Histogram(0:1/200:1), M_samples .^ -1)
        plot!(
            hist,
            st = :steps,
            title = String("M") * L"_{0\nu}^{-1} (gaussian)",
            label = "NME centred in $nme",
            legend = :topright,
            linewidth = 1.5,
            tight = true,
            fillalpha = 0.5,
        )
        xaxis!(String("M") * L"_{0\nu}^{-1}")
        yaxis!("Proability Density")
        xlims!(0, 1)
        ylims!(0, ylims()[2])

        savefig(p, "temporary_pdf.pdf")
        append_pdf!(joinpath(output, filename), "temporary_pdf.pdf", cleanup = true)
    end

    for (idx, nme) in enumerate([1.0, 1.5, 2.0, nme_belley.central, 4.0, 7.5])
        M_dist = Truncated(Normal(nme, (nme_belley.up + nme_belley.low) / 2), 0, 30)
        M_samples = rand(M_dist, length(S_samples))

        p = plot()
        hist = append!(Histogram(0:1/200:1), M_samples .^ -2)
        plot!(
            hist,
            st = :steps,
            title = String("M") * L"_{0\nu}^{-2} (gaussian)",
            label = "NME centred in $nme",
            legend = :topright,
            linewidth = 1.5,
            tight = true,
            fillalpha = 0.5,
        )
        xaxis!(String("M") * L"_{0\nu}^{-2}")
        yaxis!("Proability Density")
        xlims!(0, 1)
        ylims!(0, ylims()[2])

        savefig(p, "temporary_pdf.pdf")
        append_pdf!(joinpath(output, filename), "temporary_pdf.pdf", cleanup = true)
    end


    # m^2 pdf for gaussian sampled NMEs
    p = plot()
    mycol = "orange"

    for (idx, nme) in
        enumerate([0.5, 1.0, 1.3, 1.7, 2.0, nme_belley.central, 3.0, 4.0, 5.0, 7.5])
        M_dist = Truncated(Normal(nme, (nme_belley.up + nme_belley.low) / 2), 0, 30)
        M_samples = rand(M_dist, length(S_samples))
        marginalised_mbb2_samples_Belley = mbb2(S_samples, G = phase_space, M = M_samples)

        hist = append!(Histogram(0:0.06/200:0.06), marginalised_mbb2_samples_Belley)
        plot!(
            hist,
            st = :steps,
            title = String("Posterior distribution of ") *
                    L"m_{\beta\beta}^2" *
                    String(" varying NME (gaussian sampling)"),
            label = "NME = $nme",
            linewidth = 1.5,
            tight = true,
        )
        if (quantiles)
            vline!(
                [quantile(marginalised_mbb2_samples_Belley, 0.9)],
                linestyle = :dash,
                color = tol_mk15[idx],
            )
        end

    end

    xaxis!(L"m_{\beta\beta}^2" * String(" [eV]"))
    yaxis!("Proability Density")
    xlims!(0, 0.06)
    ylims!(0, ylims()[2])

    savefig(p, "temporary_pdf.pdf")
    append_pdf!(joinpath(output, filename), "temporary_pdf.pdf", cleanup = true)


    # m^2 pdf for fixed NMEs
    p = plot()
    for (idx, nme) in
        enumerate([0.5, 1.0, 1.3, 1.7, 2.0, nme_belley.central, 3.0, 4.0, 5.0, 7.5])
        marginalised_mbb2_samples_Belley = mbb2(S_samples, G = phase_space, M = nme)

        hist = append!(Histogram(0:0.06/200:0.06), marginalised_mbb2_samples_Belley)
        plot!(
            hist,
            st = :steps,
            title = String("Posterior distribution of ") *
                    L"m_{\beta\beta}^2" *
                    String(" varying NME (fixed)"),
            label = "NME = $nme",
            #color = tol_mk15[idx],
            linewidth = 1.5,
            tight = true,
        )
        if (quantiles)
            vline!(
                [quantile(marginalised_mbb2_samples_Belley, 0.9)],
                linestyle = :dash,
                color = tol_mk15[idx],
            )
        end

    end

    xaxis!(L"m_{\beta\beta}^2" * String(" [eV]"))
    yaxis!("Proability Density")
    xlims!(0, 0.06)
    ylims!(0, ylims()[2])

    savefig(p, "temporary_pdf.pdf")
    append_pdf!(joinpath(output, filename), "temporary_pdf.pdf", cleanup = true)


    # comparison between fixed and gaussian sampled NMEs
    p = plot()
    vline!([2.6], linestyle = :dash, color = "black", label = "Belley's central value")

    gerda_T12 = 1.83
    gerda_low_limit = mbb(1 / gerda_T12 * 10, G = phase_space, M = nme_adams_up) * 1000
    gerda_upp_limit = mbb(1 / gerda_T12 * 10, G = phase_space, M = nme_adams_low) * 1000
    hspan!(
        [gerda_low_limit, gerda_upp_limit],
        fillrange = 1,
        fillalpha = 0.2,
        color = :green,
    )

    for (idx, nme) in
        enumerate([0.5, 1.0, 1.3, 1.7, 2.0, nme_belley.central, 3.0, 4.0, 5.0, 7.5])
        mycol = "orange"
        mycol2 = "navy"
        fixed_limit = fixed_limits[idx]
        limit = limits[idx]
        if idx != 1
            scatter!(
                p,
                [nme],
                [fixed_limit * 1000],
                title = L"m_{\beta\beta}" * String(" .vs. NME"),
                seriestype = :scatter,
                legend = false,
                mc = mycol,
                ms = 4,
                ma = 0.6,
            )
        else
            scatter!(
                p,
                [nme],
                [fixed_limit * 1000],
                title = L"m_{\beta\beta}" * String(" .vs. NME"),
                seriestype = :scatter,
                label = "Fixed NME",
                legend = false,
                mc = mycol,
                ms = 4,
                ma = 0.6,
            )
        end

        if idx != 1
            scatter!(
                p,
                [nme],
                [limit * 1000],
                title = L"m_{\beta\beta}" * String(" .vs. NME"),
                seriestype = :scatter,
                markershape = :diamond,
                legend = false,
                mc = mycol2,
                ms = 4,
                ma = 0.6,
            )
        else
            scatter!(
                p,
                [nme],
                [limit * 1000],
                title = L"m_{\beta\beta}" * String(" .vs. NME"),
                seriestype = :scatter,
                markershape = :diamond,
                label = "Sampled NME",
                legend = false,
                mc = mycol2,
                ms = 4,
                ma = 0.6,
            )
        end
    end

    label_fixed = "Fixed NME"
    label_sampled = "Sampled NME"
    yaxis!(L"m_{\beta\beta}" * String(" upper limit, 90% CI [meV]"))
    xaxis!("NME value")

    savefig(p, "temporary_pdf.pdf")
    append_pdf!(joinpath(output, filename), "temporary_pdf.pdf", cleanup = true)


    # comparison between fixed at 68% and gaussian sampled NMEs
    p = plot()
    vline!([2.6], linestyle = :dash, color = "black", label = "Belley's central value")

    gerda_T12 = 1.83
    gerda_low_limit = mbb(1 / gerda_T12 * 10, G = phase_space, M = nme_adams_up) * 1000
    gerda_upp_limit = mbb(1 / gerda_T12 * 10, G = phase_space, M = nme_adams_low) * 1000
    hspan!(
        [gerda_low_limit, gerda_upp_limit],
        fillrange = 1,
        fillalpha = 0.2,
        color = :green,
    )

    for (idx, nme) in
        enumerate([0.5, 1.0, 1.3, 1.7, 2.0, nme_belley.central, 3.0, 4.0, 5.0, 7.5])
        mycol = "orange"
        mycol2 = "navy"
        mycol3 = "orangered"
        fixed_limit = fixed_limits[idx] * 1000
        limit = limits[idx] * 1000
        upp_68_limit = fixed_68_upp_limits[idx] * 1000
        low_68_limit = fixed_68_low_limits[idx] * 1000
        scatter!(
            p,
            [nme],
            [limit],
            title = L"m_{\beta\beta}" * String(" .vs. NME"),
            seriestype = :scatter,
            markershape = :diamond,
            legend = false,
            mc = mycol2,
            ms = 4,
            ma = 0.6,
        )

        if low_68_limit > 0
            scatter!(
                p,
                [nme],
                [low_68_limit],
                title = L"m_{\beta\beta}" * String(" .vs. NME"),
                seriestype = :scatter,
                markershape = :dtriangle,
                legend = false,
                mc = mycol3,
                ms = 4,
                ma = 0.6,
            )
        end
        if upp_68_limit > 0
            scatter!(
                p,
                [nme],
                [upp_68_limit],
                title = L"m_{\beta\beta}" * String(" .vs. NME"),
                seriestype = :scatter,
                markershape = :utriangle,
                legend = false,
                mc = mycol3,
                ms = 4,
                ma = 0.6,
            )
        end
    end

    yaxis!(L"m_{\beta\beta}" * String(" upper limit, 90% CI [meV]"))
    xaxis!("NME value")

    savefig(p, "temporary_pdf.pdf")
    append_pdf!(joinpath(output, filename), "temporary_pdf.pdf", cleanup = true)


    # comparison between fixed at 68% and gaussian sampled NMEs and fixed NMEs to central value
    p = plot()
    vline!([2.6], linestyle = :dash, color = "black", label = "Belley's central value")

    gerda_T12 = 1.83
    gerda_low_limit = mbb(1 / gerda_T12 * 10, G = phase_space, M = nme_adams_up) * 1000
    gerda_upp_limit = mbb(1 / gerda_T12 * 10, G = phase_space, M = nme_adams_low) * 1000
    hspan!(
        [gerda_low_limit, gerda_upp_limit],
        fillrange = 1,
        fillalpha = 0.2,
        color = :green,
    )

    for (idx, nme) in
        enumerate([0.5, 1.0, 1.3, 1.7, 2.0, nme_belley.central, 3.0, 4.0, 5.0, 7.5])
        nme_low = (nme - (nme_belley.up + nme_belley.low) / 2 / 10)
        nme_up = (nme + (nme_belley.up + nme_belley.low) / 2 / 10)

        mycol = "orange"
        mycol2 = "navy"
        mycol3 = "orangered"
        fixed_limit = fixed_limits[idx] * 1000
        limit = limits[idx] * 1000
        upp_68_limit = fixed_68_upp_limits[idx] * 1000
        low_68_limit = fixed_68_low_limits[idx] * 1000

        if low_68_limit > 0 && upp_68_limit > 0
            plot!(
                [nme_low, nme_up, nme_up, nme_low, nme_low],
                [upp_68_limit, upp_68_limit, low_68_limit, low_68_limit, upp_68_limit],
                fillalpha = 0.3,
                label = "",
                fillrange = 0,
                linewidth = 0,
                color = mycol3,
            )
        else
            if low_68_limit > 0
                plot!([nme_low, nme_up], [low_68_limit, low_68_limit], color = mycol3)
            end
            if upp_68_limit > 0
                plot!([nme_low, nme_up], [upp_68_limit, upp_68_limit], color = mycol3)
            end
        end


        scatter!(
            p,
            [nme],
            [limit],
            title = L"m_{\beta\beta}" * String(" .vs. NME"),
            seriestype = :scatter,
            markershape = :diamond,
            legend = false,
            mc = mycol2,
            ms = 4,
            ma = 0.6,
        )

        scatter!(
            p,
            [nme],
            [fixed_limit],
            title = L"m_{\beta\beta}" * String(" .vs. NME"),
            seriestype = :scatter,
            markershape = :diamond,
            legend = false,
            mc = mycol,
            ms = 4,
            ma = 0.6,
        )


    end

    yaxis!(L"m_{\beta\beta}" * String(" upper limit, 90% CI [meV]"))
    xaxis!("NME value")

    savefig(p, "temporary_pdf.pdf")
    append_pdf!(joinpath(output, filename), "temporary_pdf.pdf", cleanup = true)



    # m pdf for gaussian sampled NMEs, changing width
    p = plot()
    limits = []
    mycol = "orange"
    sigma_limits = []

    sigma_og = (nme_belley.up + nme_belley.low) / 2
    values = []
    for (idx, k) in enumerate([
        -sigma_og,
        -sigma_og / 1.5,
        -sigma_og / 2,
        -sigma_og / 2.5,
        -sigma_og / 5,
        -sigma_og / 10,
        0,
        sigma_og / 10,
        sigma_og / 5,
        sigma_og / 2.5,
        sigma_og / 2,
        sigma_og / 1.5,
        sigma_og,
    ])
        nme_sigma = sigma_og - k
        M_dist = Truncated(Normal(nme_belley.central, nme_sigma), 0, 30)
        M_samples = rand(M_dist, length(S_samples))
        marginalised_mbb_samples_Belley = mbb(S_samples, G = phase_space, M = M_samples)
        append!(sigma_limits, quantile(marginalised_mbb_samples_Belley, 0.9))
        append!(values, nme_sigma)

        sigma_limit = sigma_limits[idx]
        scatter!(
            p,
            [nme_sigma],
            [sigma_limit * 1000],
            title = L"m_{\beta\beta}" * String(" .vs. NME"),
            seriestype = :scatter,
            legend = false,
            mc = "navy",
            ms = 4,
            ma = 0.6,
        )
    end

    yaxis!(L"m_{\beta\beta}" * String(" upper limit, 90% CI [meV]"))
    xaxis!("NME sigma (NME central value: 2.60)")

    savefig(p, "temporary_pdf.pdf")
    append_pdf!(joinpath(output, filename), "temporary_pdf.pdf", cleanup = true)


end

function test_Belleys_numbers(file_path::String)

    samples = bat_read(file_path).result
    S_samples = get_S_posterior(samples)

    l1000_samples = Exponential((1 / 13) / 2.3)
    l1000_samples = rand(l1000_samples, length(S_samples))

    M_dist =
        Truncated(Normal(nme_belley.central, (nme_belley.up + nme_belley.low) / 2), 0, 30)
    M_samples = rand(M_dist, length(S_samples))

    fixed_belley_mbb_samples = mbb(1 / (13), G = phase_space, M = nme_belley.central)
    @info "L100 S,M fixed ", fixed_belley_mbb_samples * 1000

    fixed_belley_mbb_samples = mbb(l1000_samples, G = phase_space, M = M_samples)
    @info "L100 S,M sampled ", quantile(fixed_belley_mbb_samples, 0.9) * 1000

    fixed_belley_mbb_samples = mbb(1 / (13), G = phase_space, M = M_samples)
    @info "L100 S fixed, M sampled ", quantile(fixed_belley_mbb_samples, 0.9) * 1000

    fixed_belley_mbb_samples = mbb(l1000_samples, G = phase_space, M = nme_belley.central)
    @info "L100 S sampled, M fixed ", quantile(fixed_belley_mbb_samples, 0.9) * 1000

    gerda_samples = Exponential(1 / 1.83 * 10 / 2.3)
    gerda_samples = rand(gerda_samples, length(S_samples))

    fixed_belley_mbb_samples = mbb(1 / 1.83 * 10, G = phase_space, M = nme_belley.central)
    @info "GERDA S,M fixed ", fixed_belley_mbb_samples * 1000

    fixed_belley_mbb_samples = mbb(S_samples, G = phase_space, M = M_samples)
    @info "GERDA S,M sampled ", quantile(fixed_belley_mbb_samples, 0.9) * 1000

    fixed_belley_mbb_samples = mbb(1 / 1.83 * 10, G = phase_space, M = M_samples)
    @info "GERDA S fixed, M sampled ", quantile(fixed_belley_mbb_samples, 0.9) * 1000

    fixed_belley_mbb_samples = mbb(S_samples, G = phase_space, M = nme_belley.central)
    @info "GERDA S sampled, M fixed ", quantile(fixed_belley_mbb_samples, 0.9) * 1000

end
