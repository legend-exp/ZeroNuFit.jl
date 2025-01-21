using Pkg
Pkg.activate(".") 
import Pkg; Pkg.add("ArgParse")
using ArgParse
using JSON
using DataStructures
using PropDicts
using BAT
import HDF5
using TypedTables
using Plots, LaTeXStrings
using Printf
using ColorSchemes
using Random, LinearAlgebra, Statistics, Distributions, StatsBase

# style
default(
    framestyle = :box,               # Grid line transparency
    background_color = :white,       # Background color of the plot,
    titlefontsize = 10,     # Global title font size
    guidefontsize = 10,     # Global axis label font size
    tickfontsize = 10,      # Global tick label font size
    legendfontsize = 10,     # Global legend font size
    linewidth=1.5,
)
tol_colors = ColorSchemes.tol_muted

# inputs
nme = (central = 2.60, up = 1.28, low = 1.36)
phase_space = 0.237*10^-14
gA = 1.27
file_path ="../outputs/samples.h5"

# useful functions

function mbb(S; G = 1, M = 1, gA=1.27, me=511*1000)
    return sqrt.(S).*M.^-1*sqrt(10.0 ^ -27 * me^2  * G ^ -1 * gA ^ -4)
end


function mbb2(S; G = 1, M = 1, gA=1.27, me=511*1000)
    return S.*M.^-2*(10.0 ^ -27 * me^2  * G ^ -1 * gA ^ -4)
end

function get_S_posterior(samples)
    S_samples= []

    for samp in samples
        v = samp.v
        weight = samp.weight
        for w = 1:1:weight
            append!(S_samples, v[:S])
        end
    end
    S_samples
end

function plot_posterior(name,post,range,axis_label,label,quantiles;line =nothing, post_other = nothing,label_other=nothing)
    
    # plot the first histogram
    hist = append!(Histogram(0:range/200:range),post)
    p=plot(hist,
        st = :steps, 
        title = "Posterior distribution of $name",label=label,color=tol_colors[1],linewidth=1.5
    )
    if (quantiles)
        vline!([quantile(post, 0.9)],linestyle=:dash,color=tol_colors[1],label="")
    end
    if (line!=nothing)
        vline!([line],linestyle=:dash,color=tol_colors[4],label="Fixed S")
    end
    # plot the second histogram
    if (post_other!=nothing)
        hist_other = append!(Histogram(0:range/200:range),post_other)

        plot!(hist_other,
            st = :steps, 
            title = "Posterior distribution",label=label_other,color="orange"
        )
        if (quantiles)
            vline!([quantile(post_other, 0.9)],linestyle=:dash,color="orange",label="")
        end
    end
    xaxis!(axis_label)
    yaxis!("Proability Density")
    xlims!(0,range)
    ylims!(0, ylims()[2])
    display(p)

end


# get the relevant samples
samples = bat_read(file_path).result
S_samples = get_S_posterior(samples)

M_dist = Truncated(Normal(nme.central,nme.up),0,10)
M_samples = rand(M_dist, length(S_samples))

# fixed NME
fixed_mbb_samples = mbb(S_samples,G=phase_space, M = nme.central)
fixed_mbb2_samples = mbb(S_samples,G=phase_space, M = nme.central)

# vary the NME
marginalised_mbb_samples = mbb(S_samples,G=phase_space, M = M_samples)
marginalised_mbb2_samples = mbb2(S_samples,G=phase_space, M = M_samples)

# fix S and vary NME
fix_S_samples =fill(20, 60000)
M_samples = rand(M_dist, length(fix_S_samples))
fixed_S_marginalised_mbb_samples = mbb(fix_S_samples,G=phase_space, M = M_samples)
fixed_S_marginalised_mbb2_samples = mbb2(fix_S_samples,G=phase_space, M = M_samples)

fixed_S_mbb_samples = mbb(fix_S_samples,G=phase_space, M = nme.central)[1]
fixed_S_mbb2_samples = mbb2(fix_S_samples,G=phase_space, M = nme.central)[1]


# plot Mov
plot_posterior("S",S_samples,20,string("S [") * L"10^{-27} yr^{-1}]","",true)
display()
plot_posterior(L"\sqrt{S}",S_samples.^0.5,10, L"\sqrt{S} 10^{-27} yr^{-1/2}","",true)

plot_posterior("M",M_samples,10,string("M") * L"_{0\nu}","",false,line=nme.central)
plot_posterior(L"M^{-1}",M_samples.^-1,1,string("M") * L"_{0\nu}^{-1}","",false,line=1/nme.central)
plot_posterior(L"M^{-2}",M_samples.^-2,1,string("M") * L"_{0\nu}^{-2}","",false,line=1/(nme.central^2))

# plot mbb
plot_posterior(L"m_{\beta\beta}",fixed_mbb_samples,1.5,string("m") * L"_{\beta\beta}" * string(" [eV]"),"Fixed NME",true,post_other =marginalised_mbb_samples,label_other="Marginalised NME" )

# mbb^2
plot_posterior(L"m_{\beta\beta}",fixed_mbb_samples,1.5,string("m") * L"_{\beta\beta}^2" * string(" [eV")*L"^2]","Fixed NME",true,line = fixed_S_mbb_samples,post_other =marginalised_mbb_samples,label_other="Marginalised NME" )

# plot mbb
plot_posterior(L"m_{\beta\beta}",fixed_S_marginalised_mbb_samples,1.5,string("m") * L"_{\beta\beta}" * string(" [eV]"),"Marginalised NME",true,line = fixed_S_mbb_samples,post_other =nothing,label_other=nothing )

# mbb2
plot_posterior(L"m_{\beta\beta}^2",fixed_S_marginalised_mbb2_samples,.7,string("m") * L"_{\beta\beta}^2" * string(" [eV") *L"^2]","Marginalised NME",true,line = fixed_S_mbb2_samples,post_other =nothing,label_other=nothing )


sleep(100000)