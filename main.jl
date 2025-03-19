### main.jl
### -> gets a config.json file in input for running the Bayesian unbinned fit
###


using Pkg
Pkg.activate(".") # activate the environment
Pkg.instantiate() # instantiate the environment
using ArgParse
using Logging, LoggingExtras
using JSON
using FilePathsBase
# load the script to run the analysis
include("src/ZeroNuFit.jl")
using .ZeroNuFit



# process parsed arguments for the main function
function get_argparse()
    """
    Parse the script arguments
    """
    settings = ArgParseSettings(
        prog = "LEGEND ovbb Bayesian unbinned fit",
        description = "",
        commands_are_required = true,
    )
    @add_arg_table! settings begin
        "--config", "-c"
        help = "path to config file"
        arg_type = String
        required = true
    end

    parse_args(settings)
    return parse_args(settings)
end

function main()

    # read parsed arguments
    @info "running using ", Base.Threads.nthreads(), " threads"
    parsed_args = get_argparse()

    # read config path
    config_path = parsed_args["config"]
    @info "Reading configuration from: $config_path"
    config = ZeroNuFit.Utils.read_config(config_path)

    # load the output path and create the neccesary
    output_path = config["output_path"]

    for dir in [
        "$output_path/",
        "$output_path/plots/",
        "$output_path/mcmc_files/",
        "$output_path/logs/",
    ]
        if !isdir(dir)
            mkpath(dir)
        end
    end

    ZeroNuFit.Utils.set_logger(config, output_path)

    # Call the analysis function from ZeroNuFit
    ZeroNuFit.Analysis.run_analysis(config, output_path = output_path, toy_idx = nothing)

end

# Run the main function if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
