# Tutorial

The aim of this tutorial consists in building proper config JSON files in order to run a neutrinoless double-beta decay analysis over GERDA and MAJORANA DEMONSTRATO (MJD) published data.
Additional info on the meaning of input parameters can be found under the "Configuration file" section, and for input files under the "Partitions and events" section.

Table of contents:

```@contents
Pages = ["tutorial.md"]
Depth = 3
```

## GERDA Phase I: S+B fit
Let's start by fitting data acquired by the GERDA experiment during its Phase I and let's fit them with a signal+background model.

### Input files
First of all, we have to populate partitions and events JSON input files.
These dictionaries can be built manually by each user, but some reference files are already present for GERDA/LEGEND/MJD experiments at a [private location](https://github.com/tdixon97/legend-0vbb-config/tree/master) (ask [Toby Dixon](https://github.com/tdixon97) for access to the repository).
Some of these files were also copied in the [`inputs/` folder](https://github.com/legend-exp/ZeroNuFit.jl/tree/main/inputs) for already published material.
The following files can be used for this tutorial:

- GERDA Phase I Events -> `"inputs/events_gerda_pI.json"`
- GERDA Phase I Partitions -> `"inputs/partitions_gerda_pI.json"`

Notice that in `"inputs/partitions_gerda_pI.json"` we already specify in which range we want to fit data (here it is set common to all fit groups, i.e. [1930,2099] U [2109,2114] U [2124,2190] keV) and how we want to group different detectors.
In particular, 4 fit groups are specified, with names `ph1_golden`, `ph1_silver`, `phI_bege` and `phI_extra`.
For each group, we will have a separate background index (BI): $\mathcal{B}_{phI-golden}$, $\mathcal{B}_{phI-silver}$, $\mathcal{B}_{phI-bege}$, $\mathcal{B}_{phI-extra}$.
In case you want one BI only for all phase I events, $\mathcal{B}_{phI-all}$, then you have to modify the partitions input file to account for that.

### Fit configuration
Let's fit GERDA Phase I events with a simple signal+background=S+B model (`"bkg_only": false`) where the signal is modelled with a Gaussian function (default option).

The signal prior is taken as uniform (`"signal": {"prior": "uniform", ...}`) in $[0;10^{-24}]$ yr$^{-1}$.
Notice that for the signal the values are expressed in terms of $10^{-27}$yr$^{-1}$ (that is why `"signal": {"upper_bound": 1000, ...}`).
The BI prior is taken as uniform (`"bkg": {"prior": "uniform", ...}`) in $[0;0.1]$ counts/keV/kg/yr. 
Signal and background are always defined as positive, i.e. the lower bound is set to 0 by default.

We assume the 4 BIs are not correlated (`"bkg": {"correlated": {"mode": "none", "range": "none"}}`).
If they are, change the entry into `"bkg": {"correlated": {"mode": x, "range": [...,...]}}` where `x={"lognormal", "normal"}`.
These are the only _hierarchical models_ present at the moment.
More documentation on this topic can be found in "A. Gelman, J.B. Carlin, H.S. Stern, D.B. Dunson, A. Vehtari and D.B. Rubin, “Bayesian Data Analysis (3rd ed.)”, Chapman & Hall (2013)".

As regards the nuisance parameters, we can leave free the energy biases and energy widths by not fixing them to their best value (`"nuisance": {"energy_scale": {"fixed": false, ...}, ...}`), but constraining them via a Gaussian prior.
An additional option for treating all energy biases or widths together via one common parameter, i.e. $\alpha_{\Delta}$ or $\alpha_{\omega}$, can be enabled/disabled. For the moment, we leave this out and we treat nuisance parameters individually (`"nuisance": {"energy_scale": {"correlated": false, ...}, ...}`).

As regards the efficiencies, we don't fix the values to their best value (`"nuisance": {"efficiency": {"fixed": false, ...}, ...}`), but we correlated them via a global parameter $\alpha_{\varepsilon}$ (`"nuisance": {"efficiency": {"correlated": true, ...}, ...}`).

All these fit settings can therefore be grouped in the following config JSON file with name `config_gerda_phI.json`: 

```json
{
    "debug":false,
    "events":    ["inputs/events_gerda_pI.json"],
    "partitions":["inputs/partitions_gerda_pI.json"],
    "output_path": "output/fit_gerda_phI/",
    "overwrite": true,
    "light_output": false,
    "bat_fit": {"nsteps": 1e6, "nchains": 6},
    "plot": {"fit_and_data": false, "bandfit_and_data": false, "scheme":"green", "alpha": 0.3},
    "bkg_only": false,
    "signal": {"upper_bound": 1000, "prior": "uniform"},
    "bkg": {"units": "ckky", "upper_bound": 0.1, "prior": "uniform", "correlated": {"mode": "none", "range": "none"}},
    "nuisance": { 
        "energy_scale" : {
            "correlated": false,
            "fixed": false
        },
        "efficiency" : {
            "correlated": true,
            "fixed": false
        }
    }
}
```

You can now run this fit by running
```
$ julia main.jl -c config_gerda_phI.json
```


## GERDA Phase I and II: S+B fit
Let's add data acquired by the GERDA experiment during its Phase II and let's fit them with a S+B model, leaving previous settings invariate.

### Input files
In the partitions file, we group detectors from Phase II in one fit group only (i.e. `all_phase_II`) such that they all share one BI, i.e. $\mathcal{B}_{phII-all}$.

- GERDA Phase II events -> `"inputs/events_gerda_pII.json"`
- GERDA Phase II partitions -> `"inputs/partitions_gerda_pII.json"`

### Fit configuration
Combined fits of different phases of one experiment (or, from a more general point of view, combined fits of different experiments) can be achieved by introducing new events and partitions files to the config JSON file. 
Let's create the following config JSON file with name `config_gerda_phIandphII.json`: 

```json
{
    ...,
    "events":    ["inputs/events_gerda_pI.json", "inputs/events_gerda_pII.json"],
    "partitions":["inputs/partitions_gerda_pI.json", "inputs/partitions_gerda_pII.json"],
    "output_path": "output/fit_gerda_phIandphII/",
    ...
}
```
All the rest can be left unchanged and you can now run the GERDA Phase I+II combined fit by running
```
$ julia main.jl -c config_gerda_phIandphII.json
```


## GERDA Phase I and II: S+B fit, non flat background
The background is assumed flat by default, but how can we include a potential different shape?
A linear or exponential background shape were implemented as well, and one of these shapes can be specified under the `"bkg"` key.
Let's create the following config JSON file with name `config_gerda_phIandphII_linearB.json`: 

```json
{
    ...,
    "events":    ["inputs/events_gerda_pI.json", "inputs/events_gerda_pII.json"],
    "partitions":["inputs/partitions_gerda_pI.json", "inputs/partitions_gerda_pII.json"],
    "output_path": "output/fit_gerda_phIandphII_linearB/",
    "bkg": {
        "units": "ckky", 
        "upper_bound": 0.1, 
        "prior": "uniform",
        "correlated": {"mode": "none", "range": "none"},
        "shape":{
            "name":"linear",
            "pars":{"slope":[-1,3]}
        }
    },
    ...
}
```

or `config_gerda_phIandphII_expoB.json`:

```json
{
    ...,
    "events":    ["inputs/events_gerda_pI.json", "inputs/events_gerda_pII.json"],
    "partitions":["inputs/partitions_gerda_pI.json", "inputs/partitions_gerda_pII.json"],
    "output_path": "output/fit_gerda_phIandphII_expoB/",
    "bkg": {
        "units": "ckky", 
        "upper_bound": 0.1, 
        "prior": "uniform",
        "correlated": {"mode": "none", "range": "none"},
        "shape":{
            "name":"exponential",
            "pars":{"slope":[-10,10]}
        }
    },
    ...
}
```


You can now run this fit by running
```
$ julia main.jl -c config_gerda_phIandphII_linearB.json
```
if you want to shape the background with a linear function, or 
```
$ julia main.jl -c config_gerda_phIandphII_expoB.json
```
if you want to shape the background with an exponential function.



## GERDA Phase I and II: B only fit
Sometimes it is required to fit under the $\mathcal{S}=0$ (no signal) assumption.
This might be helpful when performing sensitivity studies in the context of $0\nu\beta\beta$ decay analyses.
This can be achieved by setting `"bkg_only": true` in our config JSON file that we now call `config_gerda_phIandphII_NoSignal.json`:

```json
{
    ...,
    "events":    ["inputs/events_gerda_pI.json", "inputs/events_gerda_pII.json"],
    "partitions":["inputs/partitions_gerda_pI.json", "inputs/partitions_gerda_pII.json"],
    "output_path": "output/fit_gerda_phIandphII_NoSignal/",
    "bkg_only": true,
    ...
}
```
You can now run this fit by running
```
$ julia main.jl -c config_gerda_phIandphII_NoSignal.json
```

Additional details on the type of available sensitivity studies (e.g. how to generate fake spectra and fit them) can be found in the "Generating toys" section.

## GERDA + MJD: different S shapes
As explained in ["I. J. Arnquist et al., Final Result of the Majorana Demonstrator’s Search for Neutrinoless Double- β Decay in Ge 76, PRL 130, 062501 (2023)"](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.130.062501), MJD used a modified Gaussian signal peak shape.
The code can take care of this difference once you specify the type of signal shape one wants to use in the partitions file.
Use `"signal_name": "gaussian_plus_lowEtail"` for MJD, and `"signal_name": "gaussian"` (or don't enter any key - this is the default option) for GERDA partitions.
Below we report an example of how the MJD partitions JSON file should look like for changing the Gaussian signal function into a modified one.
Notice that for MJD we can also specify a different fit range (that now includes an additional window in 2209.1-2350.0 keV) compared to the one used by GERDA.

```json
    {
        "fit_groups": {
            "mjd-DS0": {
                "range": [[1950.0, 2098.511], [2108.511, 2113.513], [2123.513, 2199.1], [2209.1, 2350.0]],
                "model": "uniform",
                "bkg_name": "B_mjd-DS0",
                "signal_name": "gaussian_plus_lowEtail"
            }
        },
        "partitions": {
            "mjd-DS0": [...]
        }
    }
```

