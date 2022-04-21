# PALEOcopse.jl documentation

The PALEOcopse.jl repository provides:

- A set of Reactions for COPSE components (land, ocean, ...). These can be made available for use by other models by adding the PALEOcopse package, eg

    julia> ] add PALEOcopse
    julia> import PALEOcopse

- Examples for the standalone COPSE reloaded model.  To use these, clone the github repository to a local folder,

    > git clone https://github.com/PALEOtoolkit/PALEOcopse.jl.git PALEOcopse
then activate the `PALEOcopse\examples` Julia environment, setting the source for the `PALEOcopse` package as the downloaded PALEOcopse repository:

    julia> cd("PALEOcopse/examples/COPSE")
    julia> import Pkg
    julia> Pkg.activate("..") # use the PALEOcopse/examples environment
    julia> Pkg.dev("../..")   # use the local version of PALEOcopse packages to allow local modifications
    julia> Pkg.instantiate() 

PALEO documentation follows the recommendations from <https://documentation.divio.com/>

