# PALEOcopse.jl documentation

## Installation and running the model

### Installation

NB: requires Julia 1.6 or later.  To check the Julia version:

    julia> versioninfo()

Clone this github repository to local directory `PALEOcopse`: from a linux bash prompt or a Windows terminal,

    $ git clone https://github.com/PALEOtoolkit/PALEOcopse.jl.git PALEOcopse

Start julia and navigate to the `PALEOcopse/examples` folder, and run `setup.jl` to configure the `PALEOcopse/examples`
Julia environment to use the local (downloaded) version of the PALEOcopse package:

    julia> cd("PALEOcopse/examples")
    julia> include("setup.jl") # use the local version of PALEOcopse packages to allow local modifications
   
### Running the model
Start julia and navigate to the `PALEOcopse` folder, then:

    julia> cd("examples/COPSE")
    julia> import Pkg
    julia> Pkg.activate("..") # use the PALEOcopse/examples environment

    julia> include("COPSE_reloaded_reloaded.jl")

## Evaluation data

Datasets are not part of the public release but are available on request from the authors.

## Using PALEOcopse Reactions from other models

The PALEO Reactions comprising the COPSE model are available when the registered PALEOcopse package is loaded (without downloading the repository), ie

    julia> Pkg.add("PALEOcopse")
    julia> import PALEOcopse


