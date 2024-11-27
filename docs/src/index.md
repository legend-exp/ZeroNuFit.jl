# ZeroNuFit.jl Documentation

Welcome to the documentation for `ZeroNuFit.jl`.

## Introduction

ZeroNuFit.jl is a Julia package for running an unbinned fit of a Gaussian signal over a background for the neutrinoless double-beta decay ($0\nu\beta\beta$) analysis.
The tool was developed for the LEGEND experiment but it can be easily extended to data collected by any other $0\nu\beta\beta$ experiment or to any other physical processes that can be modelled in a similar way.

The package uses the [BAT.jl](https://bat.github.io/BAT.jl/dev/) tool box  for Bayesian inference (see _O. Schulz, F. Beaujean, A. Caldwell, C. Grunwald, V. Hafych, K. Kröninger et al., “Bat.jl: A julia-based tool for bayesian inference”, SN Computer Science 2 (2021) 210_).

### Likelihood Function
The implemented unbinned Likelihood function reads as: $\mathcal{L} = \prod_{\rm k=1}^{N_{\rm p}} \left[\frac{\mu_{\rm k}^{N_{\rm k}} e^{-\mu_{\rm k}}}{N_{\rm k}!} \times \prod_{\rm i=1}^{N_{\rm k}} \frac{1}{\mu_{\rm k}}\left(\frac{\mu_{\rm b,\,k}}{\Delta E} + \mu_{\rm s,\,k} \times \frac{dP(E_{\rm i} | Q_{\beta\beta} +\Delta_{\rm k}, \omega_{\rm k})}{dE}\right) \right]$.

Here, the first product runs over the number of partitions _k_ ($N_{\rm p}$ partitions in total) and the second over the events _i_ in a given partition ($N_{\rm k}$ events in total).
In particular, $E_{\rm i}$ are the energies of the events falling in the analysis window for a partition _k_ with energy bias $\Delta_{\rm k}$ and energy width $\omega_{\rm k}$.

The background and signal contributions are defined as $\mu_{\rm b,\,k} = \mathcal{B}_{\rm k} \cdot \Delta E \cdot \mathcal{E}_{\rm k}$ and $\mu_{\rm s,\,k} = \frac{\text{ln}\,2\mathcal{N}_{\rm A}}{m_{\rm 76}} \cdot (\varepsilon_{\rm k} + \alpha \cdot \sigma_{\varepsilon_{\rm k}}) \cdot \mathcal{E}_{\rm k} \cdot \mathcal{S}$.
The total contribution can thus be defined as $\mu_{\rm b,\,k} = \mu_{\rm b,\,k} + \mu_{\rm s,\,k}$.

The prior distributions were taken as Gaussian distributions centred around the true value $\pi(\Delta_{\rm k},\, \omega_{\rm k},\, \alpha) = \frac{1}{\sqrt{2\pi}\sigma_{\Delta_{\rm k}}} \,e^{-\frac{\left(\Delta_{\rm k} - \widehat{\Delta}_{\rm k}  \right)^2}{2 \sigma_{\Delta_{\rm k}} ^2}} \times \frac{1}{\sqrt{2\pi}\sigma_{\omega_{\rm k}}} \,e^{-\frac{\left(\omega_{\rm k} - \widehat{\omega}_{\rm k}  \right)^2}{2 \sigma_{\omega_{\rm k}} ^2}} \times e^{-\alpha^2/2}$.

The signal energy distribution for each partition is taken as a Gaussian (e.g. for GERDA/LEGEND), $\frac{dP(E | x, \sigma)}{dE} = \frac{1}{\sqrt{2\pi\sigma^2}} \times e^{-\frac{\left(E - x \right)^2}{2 \sigma ^2}}$.

Alternatively, the signal energy distribution can also be shaped as a Gaussian with a tail at low energies (e.g. for MAJORANA DEMONSTRATOR), $\frac{dP(E | x, \gamma)}{dE} = \frac{1-f}{\sqrt{2\pi(\gamma \sigma)^2}} \times e^{-\frac{\left(E - x \right)^2}{2 (\gamma \sigma) ^2}}+\frac{f}{2\gamma\tau}\times e^{ \frac{(\gamma \sigma)^2}{2(\gamma\tau)^2} + \frac{E-x}{\gamma\tau} } \times erfc \left(\frac{\sigma}{\sqrt{2}\tau}+ \frac{E-x}{\sqrt{2}\gamma\sigma} \right)$.

In particular, for a given partition _k_, $\mathcal{B}_{\rm k}$ is the background index, $\mathcal{E}_{\rm k}$ is the exposure, $\Delta E$ is the fit window, $\varepsilon_{\rm k}$ is the efficiency with uncertainty $\sigma_{\varepsilon_{\rm k}}$ and $\mathcal{S}=\Gamma_{\rm 1/2}=1/T_{\rm 1/2}$ is the signal.

## Table of contents

```@contents
Pages = [
    "installation.md",
    "config.md",
    "inputs.md",
    "toys.md"
]
Depth = 1
```