### fitting.jl
#
# Authors: Sofia Calgaro, Toby Dixon
# 
###
module Fitting
using Random, LinearAlgebra, Statistics, Distributions, StatsBase
using BAT, DensityInterface, IntervalSets
using TypedTables
using Plots
using Cuba
using SpecialFunctions



function get_evidence(data, func, prior, method)
    likelihood = build_simple_likelihood(data, func)
    posterior = PosteriorMeasure(likelihood, prior)

    return bat_integrate(posterior, method)
end


function get_qbb_posterior(fit_function, samples)
    qbb = []

    for samp in samples
        v = samp.v
        weight = samp.weight
        for w = 1:1:weight
            append!(qbb, fit_function(NamedTuple(v), constants.Qbb))
        end
    end
    return qbb
end



function get_n_posterior(samples)
    ns = []

    for samp in samples
        v = samp.v
        weight = samp.weight
        for w = 1:1:weight
            append!(ns, v.n)
        end
    end
    return ns
end


struct LogFlat <: ContinuousUnivariateDistribution
    a::Float64
    b::Float64
end

Distributions.support(::LogFlat) = (0.0, Inf)
Distributions.pdf(d::LogFlat, x) =
    (x > d.a) && (x < d.b) ? 1 / (x * (log(d.b) - log(d.a))) : 0
Distributions.rand(d::LogFlat) = exp(rand() * (log(d.b) - log(d.a)) + log(d.a))
Distributions.logpdf(d::LogFlat, x::Float64) =
    x > d.a && x < d.b ? -log(x) - log(log(d.b) - log(d.a)) : -Inf
Distributions.cdf(d::LogFlat, x::Float64) =
    x < d.a ? 0 : x < d.b ? (log(x) - log(d.a)) / (log(d.b) - log(d.a)) : 1
Distributions.params(d::LogFlat) = (d.a, d.b)
end