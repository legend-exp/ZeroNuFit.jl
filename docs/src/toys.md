Table of contents:

```@contents
Pages = ["toys.md"]
Depth = 3
```

## Running toys for posterior predictive distribution studies

In order to check for any mis-modeling in the fits and to further understand the derived half-life limits, we performed "sensitivity" studies. 
In a Bayesian context, the sensitivity is related to the concept of posterior predictive distributions.
The prior predictive distribution is the expected distribution of data coming from a future experiment identical to that performed and repeated under the same conditions.
This marginalizes the uncertainty in the model parameters, based on their prior distributions. 
Alternatively, the posterior predictive distribution weights the data by the posterior obtained from the analysis.
If the original data were modeled appropriately, then fake data generated under the given model should distribute similarly to the original data.

Considering the observed signal counts as an observable of the data, we can thus extract "sensitivity curves" that can be interpreted as the distribution of expected future limits derived from repeating the identical experiment.
This is distinct from a frequentist sensitivity since uncertainty on nuisance parameters are marginalized over.

A module `sensitivity.jl` is present for generating toys and running sensitivity studies. The script can be run as

```
$ julia sensitivity.jl -c config_fake_data.json -i N
```

where `N` is an integer number corresponding to the toy index.
The command can be run in an external bash script for looping over this index.

The input config file (`config_fake_data.json`) has the following entries:

```
{
    "path_to_fit": "output/fit_gerda_phIandphII_NoSignal/",
    "best_fit": false,
    "seed": null
}

```

where
- `"path_to_fit"` is the path to the already performed fit over real data;
- `"best_fit": true` if we want to fix the paramaters to the best fit;
- `"seed": null` if we want a random seed when generating fake data, otherwise you can fix it to an Int value.

Any information about the signal being included or not in the fit of real data, was saved and retrieved from the output JSON file with results.

Below, we show an example of bash file used for running sensitivity studies as multiple jobs on NERSC:

```bash
#!/bin/bash                                                                                                                                                 
#SBATCH -q regular                                                                                                                                       
#SBATCH --constraint=cpu                                                                                                                                    
#SBATCH -t 48:00:00
#SBATCH -J sens_test                                                                                                                                         
#SBATCH --output parallel.log                                                     
#SBATCH --error parallel.err  

module load parallel
module load julia
srun="srun -N 1"
parallel="parallel --delay 1 -j 128"

# run parallel jobs
$srun  $parallel "julia sensitivity.jl -c config_fake_data.json -i {1}" ::: {1..10000} &

wait
```


## Using already existing toys to test alternative models
Another way to run the code is present if, for instance, an user wants to use toy data generated according to one model but fit them with another model.
In this case, the path to the folder containing the already existing JSON files with toy data has to be provided together with the toy index:

```
julia sensitivity.jl -c config_fake_data.json -i N --path_to_toys path_to_your_toys
```

Below, an updated version of a bash file that can be used for retrieving multiple existing toy data and running sensitivity studies as multiple jobs on NERSC:

```bash
#!/bin/bash                                                                                                                                                 
#SBATCH -q regular                                                                                                                                       
#SBATCH --constraint=cpu                                                                                                                                    
#SBATCH -t 48:00:00
#SBATCH -J sens_test                                                                                                                                         
#SBATCH --output parallel.log                                                     
#SBATCH --error parallel.err             

# set the directory path to toys
path_to_toys="path_to_your_toys"
all_files=("$path_to_toys"/*.json)
full_paths=()
for file in "${all_files[@]}"; do
    if [[ -f "$file" ]]; then 
        full_paths+=("$file")
    fi
done
if [ ${#full_paths[@]} -eq 0 ]; then
    echo "The list of existing toy data is empty! Exit here."
    exit 1
else
    echo "You are going to run a fit over ${#full_paths[@]} number of already existing toys stored under $path_to_toys"
fi

# array to hold toy_idx
toy_indices=()

# loop over available fake JSON toys
for path in "${full_paths[@]}"; do
    base_name="${path%.json}"
    number_str="${base_name##*fake_data}"  
    toy_idx=$((number_str))  

    toy_indices+=("$toy_idx") 
done
echo "List of toy indices: ${toy_indices[*]}"

module load parallel
module load julia
srun="srun -N 1"
parallel="parallel --delay 1 -j 128"
$srun $parallel "julia sensitivity.jl -c config/toy_9_l200_1BI_new_data_same_bkg_noS.json -i {1}" ::: "${toy_indices[@]}"

wait
```