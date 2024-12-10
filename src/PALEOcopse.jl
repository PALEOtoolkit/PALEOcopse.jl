module PALEOcopse

import PrecompileTools

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


# gains from precompilation with Julia 1.10.0-beta2, PrecompileTools v1.2.0 are 48sec -> 27sec for include("COPSE_reloaded_reloaded.jl"))
if VERSION >= v"1.8.0" # negligible benefit from precompile prior to Julia 1.8.0
    @PrecompileTools.setup_workload begin
        # create Reactions and register methods to precompile this code

        # Putting some things in `setup` can reduce the size of the
        # precompile file and potentially make loading faster.
        
        configlist = [
            (joinpath(@__DIR__, "../examples/COPSE/COPSE_bergman2004_bergman2004_cfg.yaml"), "Bergman2004"),
            # omit as these now require PALEOocean for ReactionOceanNoTransport
            # (joinpath(@__DIR__, "../examples/COPSE/COPSE_reloaded_bergman2004_cfg.yaml"), "model1"),
            # (joinpath(@__DIR__, "../examples/COPSE/COPSE_reloaded_reloaded_cfg.yaml"), "model1"),
        ]
        @PrecompileTools.compile_workload begin
            # all calls in this block will be precompiled, regardless of whether
            # they belong to your package or not (on Julia 1.8 and higher)
            logger=Logging.NullLogger()
            # logger=Logging.ConsoleLogger()
           
            for (config_file, config_model) in configlist
                PB.run_model(config_file, config_model; call_do_deriv=true, logger)
            end

        end
    end
end

end
