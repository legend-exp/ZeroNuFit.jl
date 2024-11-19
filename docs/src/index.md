# ZeroNuFit.jl Documentation

Welcome to the documentation for `ZeroNuFit.jl`.

## Introduction

`ZeroNuFit.jl` is a Julia package for running an unbinned fit of a gaussian signal over a background for the neutrinoless double beta decay analysis.

The tool is based on the [BAT.jl](https://bat.github.io/BAT.jl/dev/) package for performing a Bayesian analysis based on the following Likelihood function:

\begin{equation}
    \mathcal{L} = \prod_{\rm k=1}^{N_{\rm p}} \left[\frac{\mu_{\rm k}^{N_{\rm k}} e^{-\mu_{\rm k}}}{N_{\rm k}!} \times 
          \prod_{\rm i=1}^{N_{\rm k}} \frac{1}{\mu_{\rm k}}
          \left(\frac{\mu_{\rm b,\,k}}{\Delta E} + 
          \mu_{\rm s,\,k} \times \frac{dP(E_{\rm i} | Q_{\beta\beta} +\Delta_{\rm k}, \omega_{\rm k})}{dE}\right) \right]
\end{equation}
      
where the first product runs over the number of partitions _k_ ($N_{\rm p}$ partitions in total) and the second over the events _i_ in a given partition ($N_{\rm k}$ events in total).
In particular, $E_{\rm i}$ are the energies of the events falling in the analysis window for a partition _k_ with energy bias $\Delta_{\rm k}$ and energy width $\omega_{\rm k}$.


The background and signal contributions are defined as 
\begin{equation}
    \mu_{\rm b,\,k} = \mathcal{B}_{\rm k} \cdot \Delta E \cdot \mathcal{E}_{\rm k},
\end{equation}

\begin{equation}
    \mu_{\rm b,\,k} = \frac{\text{ln}\,2\mathcal{N}_{\rm A}}{m_{\rm 76}} \cdot (\varepsilon_{\rm k} + \alpha \cdot \sigma_{\varepsilon_{\rm k}}) \cdot \mathcal{E}_{\rm k} \cdot \mathcal{S}\,,
\end{equation}
such that the total contribution is defined as
\begin{equation}
    \mu_{\rm b,\,k} = \mu_{\rm b,\,k} + \mu_{\rm s,\,k}. 
\end{equation}

The prior distributions were taken as Gaussian distributions centred around the true value
\begin{equation}
    \pi(\Delta_{\rm k},\, \omega_{\rm k},\, \alpha) = \frac{1}{\sqrt{2\pi}\sigma_{\Delta_{\rm k}}} \,
      e^{-\frac{\left(\Delta_{\rm k} - \widehat{\Delta}_{\rm k}  \right)^2}{2 \sigma_{\Delta_{\rm k}} ^2}} \times
      \frac{1}{\sqrt{2\pi}\sigma_{\omega_{\rm k}}} \,
      e^{-\frac{\left(\omega_{\rm k} - \widehat{\omega}_{\rm k}  \right)^2}{2 \sigma_{\omega_{\rm k}} ^2}} \times e^{-\alpha^2/2} \,. 
\end{equation}

The signal energy distribution for each partition is taken to be either a Gaussian (e.g. for GERDA/LEGEND),
\begin{equation}
    \frac{dP(E | x, \sigma)}{dE} = 
    \frac{1}{\sqrt{2\pi\sigma^2}} \times
      e^{-\frac{\left(E - x \right)^2}{2 \sigma ^2}},
\end{equation}

or Gaussian with a tail at low energies (e.g. for MAJORANA DEMONSTRATOR),
\begin{equation}
    \frac{dP(E | x, \gamma)}{dE} = 
    \frac{1-f}{\sqrt{2\pi(\gamma \sigma)^2}} \times
      e^{-\frac{\left(E - x \right)^2}{2 (\gamma \sigma) ^2}}+
      \frac{f}{2\gamma\tau}\times e^{ \frac{(\gamma \sigma)^2}{2(\gamma\tau)^2} + \frac{E-x}{\gamma\tau} } 
      \times erfc \left(\frac{\sigma}{\sqrt{2}\tau}+ \frac{E-x}{\sqrt{2}\gamma\sigma} \right).
\end{equation}

In particular, for a given partition _k_, $\mathcal{B}_{\rm k}$ is the background index, $\mathcal{E}_{\rm k}$ is the exposure, $\Delta E$ is the fit window, $\varepsilon_{\rm k}$ is the efficiency with uncertainty $\sigma_{\varepsilon_{\rm k}}$ and $\mathcal{S}=\Gamma_{\rm 1/2}=1/T_{\rm 1/2}$ is the signal.

The tool was developed for the LEGEND experiment but data coming from other experiments can be fitted as well with this tool.
      

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