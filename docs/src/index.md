# ZeroNuFit.jl Documentation

Welcome to the documentation for `ZeroNuFit.jl`.

## Introduction

ZeroNuFit.jl is a Julia package for running an extended unbinned fit of a Gaussian signal over a background for the neutrinoless double-beta decay ($0\nu\beta\beta$) analysis.
The tool was developed for the LEGEND experiment but it can be easily extended to data collected by any other $0\nu\beta\beta$ experiment or to any other physical processes that can be modelled in a similar way.

The package uses the [BAT.jl](https://bat.github.io/BAT.jl/dev/) tool box  for Bayesian inference (see _O. Schulz, F. Beaujean, A. Caldwell, C. Grunwald, V. Hafych, K. Kröninger et al., “Bat.jl: A julia-based tool for bayesian inference”, SN Computer Science 2 (2021) 210_).


## Table of contents

```@contents
Pages = [
    "likelihood.md",
    "installation.md",
    "config.md",
    "inputs.md",
    "toys.md",
    "tutorial.md",
    "api.md",
]
Depth = 1
```
