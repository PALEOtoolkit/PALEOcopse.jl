module PALEOcopse

import SnoopPrecompile

import Logging
import PALEOboxes as PB

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


# TODO gains here with Julia 1.8.1 are marginal (~20s)
if false # VERSION >= v"1.8.0" # negligible benefit from precompile prior to Julia 1.8.0
    @SnoopPrecompile.precompile_setup begin
        # create Reactions and register methods to precompile this code

        # Putting some things in `setup` can reduce the size of the
        # precompile file and potentially make loading faster.
        

        rdict = PB.find_all_reactions()
        reactionlist = ["ReactionSrSed", "ReactionLandWeatheringFluxes", "ReactionSrLand", "ReactionLandArea", "ReactionForce_LIPs",
            "ReactionForce_CK_Solar", "ReactionLandWeatheringRates", "ReactionModelBergman2004", 
            "ReactionCIsotopes", "ReactionMarineBiotaCOPSE", "ReactionOceanBurialCOPSE", "ReactionAtmOcean_A", "ReactionSedCrustCOPSE", "ReactionGlobalTemperatureBerner",
            "ReactionSrOceanfloor", "ReactionForce_CPlandrelbergman2004", "ReactionForce_UDWEbergman2004", "ReactionSeafloorWeathering",
            "ReactionSrMantleCrust", "ReactionLandBergman2004", "ReactionAtmOcean_O", "ReactionForce_spreadsheet", 
            "ReactionGlobalTemperatureCK1992", "ReactionLandBiota", "ReactionForce_Bbergman2004",
        ]

        @SnoopPrecompile.precompile_all_calls begin
            # all calls in this block will be precompiled, regardless of whether
            # they belong to your package or not (on Julia 1.8 and higher)

            Logging.with_logger(Logging.NullLogger()) do
                for r in reactionlist
                    PB.precompile_reaction(rdict, r)
                end

                # 
                PB.run_model(joinpath(@__DIR__, "../examples/COPSE/COPSE_bergman2004_bergman2004_cfg.yaml"), "Bergman2004")
                PB.run_model(joinpath(@__DIR__, "../examples/COPSE/COPSE_reloaded_bergman2004_cfg.yaml"), "model1")
                PB.run_model(joinpath(@__DIR__, "../examples/COPSE/COPSE_reloaded_reloaded_cfg.yaml"), "model1")
            end
        end
    end
end

end
