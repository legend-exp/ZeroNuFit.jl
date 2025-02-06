Table of contents:

```@contents
Pages = ["config.md"]
Depth = 3
```

## Building the configuration file
Before running the code, set the input config.json file with following entries:

```json
{
    "debug":false,
    "partitions":["config/partitions_gerda_new.json","config/partitions_l200.json","config/partitions_mjd_new.json"],
    "events":    ["config/events_gerda.json","config/events_l200.json","config/events_mjd_new_part.json"],
    "output_path": "output/fit_mjd_l200_gerda_v2/",
    "overwrite":true,
    "bat_fit": {"nchains": 6, "nsteps": 1e6},
    "plot": {
            "fit_and_data": false,
            "bandfit_and_data": false,
            "scheme":"red",
            "alpha":0.3
        },
    "bkg_only": false,
    "signal": {"upper_bound":1000, "prior": "uniform"},
    "bkg": {"upper_bound":0.1,
             "prior": "uniform",
             "correlated": {"mode": "none", "range": "none"}
             },
    ...
}
```

where
- `"debug": true` if you want to display debug output on terminal;
- `"partitions"`: list of partitions JSON inputs; it takes one entry per experiment;
- `"events"`: list of events JSON inputs; it takes one entry per experiment;
- `"output_path"`: path where to store outputs (logs, plots, mcmc results);
- `"overwrite": true` if you want to overwrite a previous fit with same `output_path`; if set to `false` but no fits were previously performed (ie there are no outputs to overwrite), the code will save the output of this fit;
- `"bat_fit"`: settings for the BAT fit, `"nchains"` numer of MCMC chains to run, `"nsteps"` number of steps for each MCMC chain;
- `"plot"`: settings for plotting; `"fit_and_data": true` plots fit line over data (and CI bands if `"bandfit_and_data": true`); `"scheme":"red"` and `"alpha":0.3` are used for customizing output appearances;
- `"bkg_only": true` if we fit assuming no signal (S=0), false otherwise;
- `"signal"`: select `"upper_bound"` for the prior and the `"prior"` shape (`uniform`, `sqrt`, ...);
- `"bkg"`: select `"upper_bound"` for the prior and the `"prior"` shape (`uniform`, ...) there are several optional keys with details given below, if these are not provided the fit defaults to a flat background without correlations.


Moreover, the config requires the following block for nuisance parameters, ie energy scale (=energy bias and resolution) and efficiency:
```json
    {
    ...
    "nuisance": { 
         "energy_bias": {
            "fixed": false,
            "correlated": false
         },
         "energy_res": {
            "fixed": false,
            "correlated": false
         },
         "efficiency" : {
            "correlated": true,
            "fixed": false
            }
    }
```

In particular, you can set `"correlated": true` if you want to use one variable to correlate the nuisance parameters (eg to speed up the computation times), and `"fixed": false` if you want to include a prior for nuisance parameters (otherwise these parameters they will be fixed to their partition value and not constrained).
 
If a variable is correlated (either `energy_bias` or `energy_res` or `efficiency`), the code will search for a field in the `fit_groups` block of the partitions JSON file to use a correlated variable per each fit group. 
In particular, the field has to be specified as:
- `"energy_bias_group_name": "..."`
- `"energy_res_group_name": "..."`
- `"efficiency_group_name": "..."`
 
Parameters are then added to the model called `αr_\$name` (for resolution), `αe_\$name` for efficiency and `αb_\$name` for bias.

!!! warning

    The $\alpha$ parameter names default to `_all`, if you want one different per experiment this must be explicitly specified in the fit groups entry

## Background shape and correlation
There are several options to control the background in more detail. These can be added to the "bkg" section of the config:
In particular:
 - "correlated" adds a hierarchical (correlated) background to the model, this key should have a dictionary giving details on the prior shape and ranges. The three options for the mode are: `none` (default), `lognormal` or `normal`. The range entry is associated to the range used for the uniform prior on the `\$sigma_B` parameter. For example:

```json
"correlated": {"mode":"lognormal","range":[0,1]}
```


- `"shape"` changes the shape of the background from uniform. The user should provide a dictionary giving details on the shape. For example:

```json
"shape": {
            "name":"exponential",
            "pars":{"slope":[-10,10]}
        },
```

The `"pars"` subdictionary describes the range of the priors on the parameters of the model, currently implemented shapes are `"uniform"`, `"linear"` and `"exponential"`. These names correspond to functions in `fitting.jl` and logical conditions in `get_bkg_pdf` in `likelihood.jl`.

This will add parameters `${bkg_name}_slope` or similar to the model (and then call them). This names therefore must correspond to the names in the functions in `fitting.jl`. To add a new shape simply define a new method in `fitting.jl` and a new logical condition in `get_bkg_pdf` in `likelihood.jl`.

!!! note

    If these keys are not provided the model defaults to a uniform uncorrelated background.

## Signal shape
The signal shape can be specified in the partition JSON file in order to have a different signal shape for different partitions.
An example is given by the combined fit of GERDA/LEGEND and MAJORANA DEMONSTRATOR (MJD) partitions where for the first one we want to model the signal with a simple gaussian function while we want to model data from MJD accounting for an additional tail at the low-energy side of the gaussian.
See the "Partitions and events" documentation.
