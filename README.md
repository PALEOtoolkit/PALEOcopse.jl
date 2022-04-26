# PALEOcopse.jl

[![CI](https://github.com/PALEOtoolkit/PALEOcopse.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/PALEOtoolkit/PALEOcopse.jl/actions/workflows/CI.yml)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://PALEOtoolkit.github.io/PALEOcopse.jl/dev)

Julia version of the COPSE (Carbon, Oxygen, Phosphorus, Sulphur and Evolution) biogeochemical model. 

The model predicts the coupled histories and controls on atmospheric O2, CO2 and ocean composition over Phanerozoic time, and is described in the following publications:

Bergman, N. M., Lenton, T. M., & Watson, A. J. (2004). COPSE: A new model of biogeochemical cycling over Phanerozoic time. American Journal of Science, 304(5), 397–437. https://doi.org/10.2475/ajs.304.5.397

Mills, B., Daines, S. J., & Lenton, T. M. (2014). Changing tectonic controls on the long-term carbon cycle from Mesozoic to present. Geochemistry, Geophysics, Geosystems, 15(12), 4866–4884. https://doi.org/10.1002/2014GC005530

Lenton, T. M., Dahl, T. W., Daines, S. J., Mills, B. J. W., Ozaki, K., Saltzman, M. R., & Porada, P. (2016). Earliest land plants created modern levels of atmospheric oxygen. Proceedings of the National Academy of Sciences, 113(35), 9704–9709. https://doi.org/10.1073/pnas.1604787113

Lenton, T. M., Daines, S.J., Mills, B. J. W. (2018). COPSE reloaded: An improved model of biogeochemical cycling over Phanerozoic time. Earth Science Reviews, 178, 1-28. https://10.1016/j.earscirev.2017.12.004

## Julia implementation

The Julia version of COPSE uses the [PALEOboxes](https://github.com/PALEOtoolkit/PALEOboxes.jl) coupler and [PALEOmodel](https://github.com/PALEOtoolkit/PALEOmodel.jl) solvers.

The COPSE model is componentized into land, atmosphere, ocean, oceanfloor and sedcrust Domains with explicit transfer of biogeochemical fluxes between them, allowing the PALEO Reactions for each Domain to be reused to create a hierarchy of coupled models eg including a spatially-resolved ocean component.

The implementation is mathematically equivalent to the [Matlab version](https://github.com/sjdaines/COPSE) and is tested against archived model output. 

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