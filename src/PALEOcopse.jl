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


# Negligible difference Julia 1.7.3
# function precompile_reactions()
#     rdict = PB.find_all_reactions()
    
#     reactionlist = ["ReactionSrSed", "ReactionLandWeatheringFluxes", "ReactionSrLand", "ReactionLandArea", "ReactionForce_LIPs",
#         "ReactionForce_CK_Solar", "ReactionLandWeatheringRates", "ReactionModelBergman2004", 
#         "ReactionCIsotopes", "ReactionOceanCOPSE", "ReactionAtmOcean_A", "ReactionSedCrustCOPSE", "ReactionGlobalTemperatureBerner",
#         "ReactionSrOceanfloor", "ReactionForce_CPlandrelbergman2004", "ReactionForce_UDWEbergman2004", "ReactionSeafloorWeathering",
#         "ReactionSrMantleCrust", "ReactionLandBergman2004", "ReactionAtmOcean_O", "ReactionForce_spreadsheet", 
#         "ReactionGlobalTemperatureCK1992", "ReactionLandBiota", "ReactionForce_Bbergman2004",
#     ]
#     for r in reactionlist
#         PB.precompile_reaction(rdict, r)
#     end

#     return nothing
# end

# precompile_reactions()

end
