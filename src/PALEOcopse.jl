module PALEOcopse

function moduledir()
    return dirname(@__DIR__)
end


include("COPSE/COPSE.jl")
include("Forcings/Forcings.jl")
include("global/Global.jl")
include("landsurface/Land.jl")
include("ocean/Ocean.jl")
include("oceanfloor/Oceanfloor.jl")
include("sedcrust/SedCrust.jl")
include("biogeochem/BioGeoChem.jl")

end
