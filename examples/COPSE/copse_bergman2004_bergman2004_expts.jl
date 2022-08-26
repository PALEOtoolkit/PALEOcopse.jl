

import PALEOboxes as PB
import PALEOmodel

function copse_bergman2004_bergman2004_expts(
    expts;   
    modelpars=Dict{}(),
)

    model = PB.create_model_from_config(
        joinpath(@__DIR__, "COPSE_bergman2004_bergman2004_cfg.yaml"), 
        "Bergman2004"; 
        modelpars=modelpars,
    )
        
    ###############################################
    # choose an 'expt' (a delta to the base model)
    ###############################################

    for expt in expts
        if expt == "baseline"
            # baseline configuration
            
        elseif length(expt)==2 && expt[1] == "tforce_constant"
            # Set constant forcing time
            # NB: use with 'modelpars'=Dict("tforcevar"=>"tforce_constant")        
            PB.set_variable_attribute!(model, "global", "tforce_constant", :initial_value, expt[2])

        else
            error("unknown expt ", expt)
        end
    end

    # retrieve COPSE model instance and set parameters for steady-state
    global_domain = PB.get_domain(model, "global")
    copsemodel = PB.get_reaction(global_domain, "model_Bergman2004")
    PALEOcopse.COPSE.ModelBergman2004.set_parameters_modern_steady_state(copsemodel)

    return model
end


function copse_bergman2004_bergman2004_plot(output; pager=PALEOmodel.DefaultPlotPager())
    
    # Conservation checks
    pager(
        plot(title="Total carbon",          output, "global.".*["total_C"], ylabel="mol C",),
        plot(title="Total carbon moldelta", output, "global.".*["total_C.v_moldelta"], ylabel="mol C * per mil",),
        plot(title="Total sulphur",         output, "global.".*["total_S"], ylabel="mol S",),
        plot(title="Total sulphur moldelta",output, "global.".*["total_S.v_moldelta"], ylabel="mol S * per mil",),
        plot(title="Total redox",           output, "global.".*["total_redox"], ylabel="mol O2 equiv",),
        :newpage,

        # Forcings
        plot(title="Physical forcings",     output, "global.".*["DEGASS", "UPLIFT"],  ylabel="normalized forcing"),
        plot(title="Evolutionary forcings", output, "global.".*["EVO", "W", "Bforcing", "CPland_relative"], ylabel="normalized forcing",),    

        # Outputs
        plot(title="Carbon",                output, "global.".*["total_C", "C", "G", "A"], ylabel="mol C",),    
        plot(title="Sulphur",               output, "global.".*["total_S", "S", "PYR", "GYP"], ylabel="mol S",),    
        plot(title="Carbon isotopes",       output, "global.".*["mccb_delta", "A_delta", "C_delta"], ylabel="delta 13C (per mil)",),

        :newpage # flush any partial page
    )

    return nothing
end

function copse_bergman2004_bergman2004_compare(outputs; pager=PALEOmodel.DefaultPlotPager())
    
    pager(
        plot(title="Oxygen",            outputs,    ["global.pO2PAL", "global.ANOX"], ylabel="pO2 (PAL), anoxia fraction"),
        plot(title="pCO2",              outputs,    "global.pCO2PAL", ylabel="pCO2 (PAL)"),
        plot(title="P",                 outputs,    "global.P", ylabel="P"),
        plot(title="N",                 outputs,    "global.N", ylabel="N"),
        plot(title="Sulphur isotopes",  outputs,    "global.S_delta", ylabel="delta 13C (‰)"),
        plot(title="Carbonate isotopes", outputs,   "global.mccb_delta", ylabel="delta 13C (‰)"),
    
        :newpage # flush any partial page
    )

    return nothing
end
