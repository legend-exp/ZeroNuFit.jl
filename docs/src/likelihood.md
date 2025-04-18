# Likelihood implementation

Table of contents:

```@contents
Pages = ["likelihood.md"]
Depth = 3
```

## Likelihood Function
The implemented unbinned Likelihood function reads as: 

```math
\begin{aligned}
    \mathcal{L}(\Gamma,\, \boldsymbol{BI},\,\boldsymbol{\theta}|D) = \prod_k \bigg[ \textrm{Pois}(s_k+b_k) \bigg[ \prod_{i_k=1}^{N_k} \frac{1}{s_k + b_k} \left( b_k\cdot p_{\rm b}(E) + s_{\rm k}\cdot p_{\rm s}(E) \right)  \bigg] \bigg]
\end{aligned}
```

where $\Gamma$ is the signal rate, BI is the background index, $\boldsymbol{\theta}$ are the nuisance parameters, and $D$ are the observed data.
Here, the first product runs over the number of partitions _k_ ($N_{\rm p}$ partitions in total) and the second over the events _i_ in a given partition ($N_{\rm k}$ events in total).
In case no events are found in a given partition _k_, the above Likelihood expression simplifies into

```math
\begin{aligned}
    \mathcal{L}(\Gamma,\, \boldsymbol{BI},\,\boldsymbol{\theta}|D) = \prod_k \textrm{Pois}(s_k+b_k) 
\end{aligned}
```

The terms $p_{\rm b}(E)$ and $p_{\rm s}(E)$ represent the background and signal distributions, respectively, both expressed as a function of energy $E$.
The background distribution can be selected among different options: flat, linear or exponential.
The signal distribution can be expressed in the following way

```math
\begin{aligned}
    p_{\rm s}(E) = \frac{dP(E_{\rm i} | Q_{\beta\beta} - \Delta_{\rm k}, \omega_{\rm k})}{dE} 
\end{aligned}
```

for a given energy $E_{\rm i}$ of an event falling in the analysis window for the partition _k_ with energy bias $\Delta_{\rm k}$ and energy width $\omega_{\rm k}$.
Here, we are assuming the energy biases were defined as $E_{\rm true} - E_{\rm cal}$.
Thus, we correct for the energy bias by adding the calculated bias to the calibrated event energy.

Taking $x=Q_{\beta\beta} - \Delta_{\rm k}$, the signal energy distribution for each partition can be taken as a Gaussian (e.g. for GERDA/LEGEND), 

```math
\begin{aligned}
\frac{dP(E | x, \sigma)}{dE} = \frac{1}{\sqrt{2\pi\sigma^2}} \times e^{-\frac{\left(E - x \right)^2}{2 \sigma ^2}}
\end{aligned}
```

Alternatively, the signal energy distribution can also be shaped as a Gaussian with a tail at low energies (e.g. for MAJORANA DEMONSTRATOR data), 

```math
\begin{aligned}
\frac{dP(E | x, \gamma)}{dE} = \frac{1-f}{\sqrt{2\pi(\gamma \sigma)^2}} \times e^{-\frac{\left(E - x \right)^2}{2 (\gamma \sigma) ^2}}+\frac{f}{2\gamma\tau}\times e^{ \frac{(\gamma \sigma)^2}{2(\gamma\tau)^2} + \frac{E-x}{\gamma\tau} } \times erfc \left(\frac{\sigma}{\sqrt{2}\tau}+ \frac{E-x}{\sqrt{2}\gamma\sigma} \right)
\end{aligned}
```

where $f$ is the fraction of events in the tail (taken as fixed), $\tau$ is the scale parameter of the tail. The $\gamma$ parameter correlates the uncertainties on both $\sigma$ and $\tau$ by simultaneously scaling $\sigma$ and $\tau$ with nominal value $\tilde{\gamma} = 1$ and uncertainty $\delta_\gamma$. Therefore, $\sigma$ and $\tau$ are taken as fixed parameters and all uncertainty is handled by $\gamma$.

The background and signal counts, i.e. $b_{\rm k}$ and $s_{\rm k}$ respectively, are defined as

```math
\begin{aligned}
b_{\rm k} = \mathcal{B}_{\rm b} \cdot \Delta E \cdot \mathcal{E}_{\rm k}
\end{aligned}
```

and

```math
\begin{aligned}
s_{\rm k} = \frac{\text{ln}\,2\cdot \mathcal{N}_{\rm A}}{m_{\rm 76}} \cdot (\varepsilon_{\rm k} + \alpha \cdot \sigma_{\varepsilon_{\rm k}}) \cdot \mathcal{E}_{\rm k} \cdot \Gamma
\end{aligned}
```

Notice that an extra index $b$ was introduced to account for different background indexes $\mathcal{B}_{\rm b}$ that might be shared across different partitions.
In particular, for a given partition _k_, $\mathcal{E}_{\rm k}$ is the exposure, $\Delta E$ is the net width of the fit window, $\varepsilon_{\rm k}$ is the efficiency with uncertainty $\sigma_{\varepsilon_{\rm k}}$ and $\Gamma$ is the signal rate.


## Likelihood implementation
In our framework, the above Likelihood product was evaluated taking the logarithm of it. 
We defined our "log Likelihood" ($LL$) as:

```math
\begin{aligned}
    LL(\Gamma)  = \underbrace{\sum_{\rm j}\sum_{\rm i_{\rm k}=1}^{N_{\rm k}} \left[\text{log}\left(Pois(b_{\rm j}+s_{\rm j})\right) + \text{log}\left(b_{\rm j} \cdot p_{\rm b}(E) + s_{\rm j} \cdot p_{\rm s}(E) \right) - \text{log}\left(b_{\rm j}+s_{\rm j} \right) \right]}_{N_{\rm k}\text{ events in $j$}} - \underbrace{\sum_{l} \left(b_{l}+s_{l} \right)}_{\text{0 events in $l$}}
\end{aligned}
```

The sum over all partitions $k$ was separated in a sum over partitions containing an event $i$ with energy $E_{\rm i}$ (sum with index $j$) and in a sum over partitions with no events (sum with index $l$).


## Prior terms

Different free prameters can be identified within the framework:
- signal, $\Gamma$
- background indeces, $\mathcal{B}_{\rm b}$ (with $b=1,\,...,\,N$)
- energy bias at $Q_{\beta\beta}$, $\Delta = E_{\textrm{true}} - E_{\textrm{cal}}$ (keV)
- energy resolution, $\sigma$ (keV)
- peak shape scale parameter, $\gamma$ (for the modified Gaussian peak shape only)
- signal efficiencies, $\varepsilon$

For the signal prior, either a uniform or a $1/\sqrt{\Gamma}$ prior can be used.

For the background indeces, a uniform prior can be used.

For nuisance parameters $\theta$ (energy biases, energy resolutions and efficiencies), Gaussian distributions centred around the true values were used

```math
\begin{aligned}
    \mathcal{P}(\boldsymbol{\theta}) =\mathcal{P}(\alpha) \cdot  \prod_{\rm k} \mathcal{P}_{\rm k}(\Delta)\cdot \mathcal{P}_{\rm k}(\sigma) =  e^{-\alpha^2/2} \cdot \prod_{\rm k} \frac{1}{\sqrt{2\pi}\sigma_{\Delta_{\rm k}}} \,
      e^{-\frac{\left(\Delta_{\rm k} - \widehat{\Delta}_{\rm k}  \right)^2}{2 \sigma_{\Delta_{\rm k}} ^2}} \cdot
      \frac{1}{\sqrt{2\pi}\sigma_{\sigma_{\rm k}}} \,
      e^{-\frac{\left(\sigma_{\rm k} - \widehat{\sigma}_{\rm k}  \right)^2}{2 \sigma_{\sigma_{\rm k}} ^2}}
\end{aligned}
```

Alternatively, the code allows the user to fix energy biases, energy resolutions and efficiencies at their true values, removing any prior constraint.
Another feature of the framework includes the possibility to correlate energy biases and resolutions via a unique term, $\alpha_{\rm bias}$ or $\alpha_{\rm reso}$, as it is done for the efficiencies in the above case.
Indeed, $\alpha$ is used as a scaling parameter that fully correlated the efficiency uncertainties across all partitions, with nominal value of 0 and standard deviation of 1.
Separating each partition efficiency into uncorrelated components and adding a nuisance parameter for each would require a great deal of CPU. 

Accounting for the different signal shape (i.e. modified Gaussian with a tail at low energies), we includeed a Gaussian prior term for the peak position offset ($\mu$) and the peak shape scale parameter ($\gamma$):

```math
\begin{aligned}
      \mathcal{P}_{\rm mod.}(\boldsymbol{\theta}) =\mathcal{P}(\alpha) \cdot   \prod_{\rm k} \mathcal{P}_{\rm k}(\mu)\cdot \mathcal{P}_{\rm k}(\sigma)=e^{-\alpha^2/2} \cdot\prod_{\rm k}  \frac{1}{\sqrt{2\pi}\sigma_{\mu_{\rm k}}} \,
      e^{-\frac{\left(\mu_{\rm k} - \widehat{\mu}_{\rm k}  \right)^2}{2 \sigma_{\mu_{\rm k}} ^2}} \cdot
      \frac{1}{\sqrt{2\pi}\sigma_{\gamma_{\rm k}}} \,
      e^{-\frac{\left(\gamma_{\rm k} - \widehat{\gamma}_{\rm k}  \right)^2}{2 \sigma_{\gamma_{\rm k}} ^2}}  
\end{aligned}
```

In the code, the product of priors was further simplified by introducing a prior term only for those partitions $j$ for which there is at least an event. 
The above products, then, can be expressed again as

```math
\begin{aligned}
    \mathcal{P}(\boldsymbol{\theta}) =\mathcal{P}(\alpha) \cdot  \prod_{\rm j} \mathcal{P}_{\rm j}(\Delta)\cdot \mathcal{P}_{\rm j}(\sigma) \text{\, and \,}
    \mathcal{P}_{\rm mod.}(\boldsymbol{\theta}) =\mathcal{P}(\alpha) \cdot  \prod_{\rm j} \mathcal{P}_{\rm j}(\mu)\cdot \mathcal{P}_{\rm j}(\sigma)
\end{aligned}
```

### Alternative background shapes

The default background modeling shape is a flat function:

```math
\begin{aligned}
  f_{\rm flat}(E) = b_{\rm k} \cdot p_{\rm b,\,flat}(E) = BI \cdot \Delta E \cdot \mathcal{E}_{\rm k} \cdot \underbrace{\frac{1}{K_{\rm flat}}}_{p_{\rm b,\,flat}(E)}
\end{aligned}
```

where $E$ is the energy and $K_{\rm flat}$ is a normalization factor that accounts for the net width of the fit window, accounting for any removed gap within it.

In the linear case, we can model the background as

```math
\begin{aligned}
  f_{\rm lin}(E) = b_{\rm k} \cdot p_{\rm b,\,lin}(E) = BI \cdot \Delta E \cdot \mathcal{E}_{\rm k} \cdot \underbrace{\left(1+ \frac{m_{\rm lin} \cdot (E-E_{\rm 0})}{E_{\rm max}-E_{\rm 0}} \right)\cdot \frac{1}{K_{\rm lin}}}_{p_{\rm b,\,lin}(E)}
\end{aligned}
```

where $m_{\rm lin}$ is the slope of the linear function, $E_{\rm 0}$ ($E_{\rm max}$) is the starting (ending) energy value of the fit window, and $K_{\rm lin}$ is the normalization factor.

In the exponential case, we can model the background as

```math
\begin{aligned}
  f_{\rm exp}(E) = b_{\rm k} \cdot p_{\rm b,\,exp}(E) = BI \cdot \Delta E \cdot \mathcal{E}_{\rm k} \cdot \underbrace{e ^ {\left(E-E_{\rm 0}\right)\cdot \frac{m_{\rm exp}}{\Delta E} } \cdot \frac{1}{K_{\rm exp}}}_{p_{\rm b,\,exp}(E)} 
\end{aligned}
```

with corresponding slope $m_{\rm exp}$ and normalization factor $K_{\rm exp}$.

The normalization factors can be expressed in a general form in the following way:

```math
\begin{aligned}
    \begin{cases}
      K_{\rm flat} = \Delta E\\[15pt]
      K_{\rm lin} = \Delta E \cdot \left(1 - \frac{m_{\rm lin}\cdot E_{\rm 0}}{E_{\rm max}-E_{\rm 0}} \right) + m_{\rm lin}\cdot \frac{\sum_{i} \left( E_{\rm h,\,i}^2 - E_{\rm l,\,i}^2 \right)}{2\left(E_{\rm max}-E_{\rm 0} \right )}\\[15pt]
      K_{\rm exp} = \frac{E_{\rm max}-E_{\rm 0}}{m_{\rm exp}} \cdot \left[\sum_{\rm i}e^{(E_{\rm h,i}-E_{\rm 0})\cdot \frac{m_{\rm exp}}{E_{\rm max}-E_{\rm 0}}} - \sum_{\rm i}e^{(E_{\rm l,i}-E_{\rm 0})\cdot \frac{m_{\rm exp}}{E_{\rm max}-E_{\rm 0}}} \right]
    \end{cases}  
\end{aligned}
```

where $E_{\rm l,i}$ ($E_{\rm h,i}$) is the starting (ending) energy value of the sub-windows defined within the analysis fit window, accounting for any excluded gap in the fit window. 
If no gaps are present, then $E_{\rm l,i}\equiv E_{\rm 0}$ and $E_{\rm h,i}\equiv E_{\rm max}$.
In particular, $\Delta E= \sum_{\rm i} \left(E_{\rm h,i}-E_{\rm l,i} \right)$.

Notice that $K_{\rm lin} \rightarrow K_{\rm flat}$ and $K_{\rm exp} \rightarrow K_{\rm flat}$ in the limit of $m_{\rm lin} \rightarrow 0$ and $m_{\rm exp} \rightarrow 0$, respectively. 
The free parameters in the linear and exponential case are $m_{\rm lin}$ and $m_{\rm exp}$ (other than the background indices, $BI$).
For all these parameters, we used a uniform prior as default, while prior ranges can be modified by the user via the configuration file. 

### Posterior distributions and marginalization

The combined posterior probability density function is calculated according to Bayes’ theorem as:

```math
\begin{aligned}
    \mathcal{P}(\Gamma,\, \boldsymbol{BI},\, \boldsymbol{\theta}|D ) \propto \underbrace{\mathcal{L}(\Gamma,\, \boldsymbol{BI},\,\boldsymbol{\theta}|D)}_{\text{Likelihood}} \cdot \underbrace{\mathcal{P}(\boldsymbol{\theta}) \cdot \mathcal{P}(\Gamma)\cdot \prod_{\rm b=1}^{N_{\rm b}} \mathcal{P}_{\rm b}(BI)}_{\text{prior terms}}
\end{aligned}
```

where $D$ are our data, i.e. the energy events surviving all cuts in the fit window. 
Here, we expressed the general case where a total number of $N_{\rm b}$ BIs are introduced. 

The marginalization is performed with the BAT toolkit via a Markov chain Monte Carlo (MCMC) numerical integration. The marginalization used the Metropolis-Hastings sampling algorithm implemented in BAT. 
The number of MCMC chains and the number of steps in each MCMC can be selected by the user.