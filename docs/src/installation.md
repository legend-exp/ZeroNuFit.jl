# First steps

## How to run the code
Run the following command by specifying the path to the configuration file used for settings:

```
$ julia main.jl -c config/config.json
```

## Julia Project enviroments
To run the code in a virtual enviroment you can use the following.
Start by entering the Julia command view by typing `julia` on your terminal.
Once inside, enter the packagge manager (via `]`), activate the environment in the current directory, and resolve (and install, if necessary) dependencies:
```julia
pkg> ] 
pkg> activate .
pkg> instantiate
```
Now you can run the script inside this enviroment with:

```
$ julia main.jl --project=. -c config/config.json
```

