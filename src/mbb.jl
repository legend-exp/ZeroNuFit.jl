using Pkg
Pkg.activate(".")
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
using Printf
using ColorSchemes
using Random, LinearAlgebra, Statistics, Distributions, StatsBase

include("constants.jl")

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

# inputs
nme = (central = constants.nme_central, up = constants.nme_up, low = constants.nme_low)
phase_space = constants.phase_space
gA = constants.gA

function mbb(S; G = 1, M = 1, gA = constants.gA, me = constants.me_keV * 1000)
    return sqrt.(S) .* M .^ -1 * sqrt(10.0^-27 * me^2 * G^-1 * gA^-4)
end


function mbb2(S; G = 1, M = 1, gA = constants.gA, me = constants.me_keV * 1000)
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

function plot_posterior(
    name::Union{String,LaTeXString},
    post,
    max_x,
    axis_label::Union{String,LaTeXString},
    label::String,
    quantiles::Bool,
    output_path::String;
    line = nothing,
    post_other = nothing,
    label_other = nothing,
)

    # plot the first histogram
    hist = append!(Histogram(0:max_x/200:max_x), post)
    p = plot(
        hist,
        st = :steps,
        title = "Posterior distribution of $name",
        label = label,
        color = tol_colors[1],
        linewidth = 1.5,
    )
    if (quantiles)
        vline!([quantile(post, 0.9)], linestyle = :dash, color = tol_colors[1], label = "")
    end
    if (line != nothing)
        vline!(
            [line],
            linestyle = :dash,
            color = tol_colors[4],
            label = String("Central M") * L"_{0\nu}^{-1}" * String("value"),
        )
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
                label = "",
            )
        end
    end
    xaxis!(axis_label)
    yaxis!("Proability Density")
    xlims!(0, max_x)
    ylims!(0, ylims()[2])

    savefig(p, "temporary_pdf.pdf")
    append_pdf!(joinpath(output, "mbb.pdf"), "temporary_pdf.pdf", cleanup = true)

end


function plot_mbb_plots(file_path::String, output_path::String)

    # remove previous mbb pdf if already present
    if isfile(joinpath(output_path, "mbb.pdf"))
        Filesystem.rm(joinpath(output_path, "mbb.pdf"), force = true)
    end

    # get the relevant samples
    samples = bat_read(file_path).result
    S_samples = get_S_posterior(samples)

    M_dist = Truncated(Normal(nme.central, nme.up), 0, 10)
    M_samples = rand(M_dist, length(S_samples))

    # fixed NME
    fixed_mbb_samples = mbb(S_samples, G = phase_space, M = nme.central)
    fixed_mbb2_samples = mbb(S_samples, G = phase_space, M = nme.central)

    # vary the NME
    marginalised_mbb_samples = mbb(S_samples, G = phase_space, M = M_samples)
    marginalised_mbb2_samples = mbb2(S_samples, G = phase_space, M = M_samples)

    # fix S and vary NME
    fix_S_samples = fill(20, 60000)
    M_samples = rand(M_dist, length(fix_S_samples))
    fixed_S_marginalised_mbb_samples = mbb(fix_S_samples, G = phase_space, M = M_samples)
    fixed_S_marginalised_mbb2_samples = mbb2(fix_S_samples, G = phase_space, M = M_samples)

    fixed_S_mbb_samples = mbb(fix_S_samples, G = phase_space, M = nme.central)[1]
    fixed_S_mbb2_samples = mbb2(fix_S_samples, G = phase_space, M = nme.central)[1]

    # plot Mov
    plot_posterior(
        "S",
        S_samples,
        30,
        string("S [") * L"10^{-27} yr^{-1}]",
        "",
        true,
        output_path,
    )
    plot_posterior(
        L"\sqrt{S}",
        S_samples .^ 0.5,
        10,
        L"\sqrt{S} 10^{-27} yr^{-1/2}",
        "",
        true,
        output_path,
    )

    plot_posterior(
        "M",
        M_samples,
        10,
        string("M") * L"_{0\nu}",
        "",
        false,
        output_path,
        line = nme.central,
    )
    plot_posterior(
        L"M^{-1}",
        M_samples .^ -1,
        1,
        string("M") * L"_{0\nu}^{-1}",
        "",
        false,
        output_path,
        line = 1 / nme.central,
    )
    plot_posterior(
        L"M^{-2}",
        M_samples .^ -2,
        1,
        string("M") * L"_{0\nu}^{-2}",
        "",
        false,
        output_path,
        line = 1 / (nme.central^2),
    )

    # Fixed NME - mbb
    plot_posterior(
        L"m_{\beta\beta}",
        fixed_mbb_samples,
        1.5,
        string("m") * L"_{\beta\beta}" * string(" [eV]"),
        "Fixed NME",
        true,
        output_path,
        post_other = marginalised_mbb_samples,
        label_other = "Marginalised NME",
    )

    # Fixed NME - mbb^2
    plot_posterior(
        L"m_{\beta\beta}",
        fixed_mbb_samples,
        1.5,
        string("m") * L"_{\beta\beta}^2" * string(" [eV") * L"^2]",
        "Fixed NME",
        true,
        output_path,
        line = fixed_S_mbb_samples,
        post_other = marginalised_mbb_samples,
        label_other = "Marginalised NME",
    )

    # Marginalised NME - mbb
    plot_posterior(
        L"m_{\beta\beta}",
        fixed_S_marginalised_mbb_samples,
        1.5,
        string("m") * L"_{\beta\beta}" * string(" [eV]"),
        "Marginalised NME",
        true,
        output_path,
        line = fixed_S_mbb_samples,
        post_other = nothing,
        label_other = nothing,
    )

    # Marginalised NME - mbb^2
    plot_posterior(
        L"m_{\beta\beta}^2",
        fixed_S_marginalised_mbb2_samples,
        0.7,
        string("m") * L"_{\beta\beta}^2" * string(" [eV") * L"^2]",
        "Marginalised NME",
        true,
        output_path,
        line = fixed_S_mbb2_samples,
        post_other = nothing,
        label_other = nothing,
    )

end
