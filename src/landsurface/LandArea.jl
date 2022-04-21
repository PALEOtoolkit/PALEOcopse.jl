module LandArea

import PALEOboxes as PB

"""
    ReactionLandArea

Calculate land areas and lithology.

Provides `GRAN_AREA` and `BA_AREA`.
"""
Base.@kwdef mutable struct ReactionLandArea{P} <:  PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
       
        PB.ParString("f_granitearea", "Fixed", 
            allowed_values=["Fixed", "Forced", "OrgEvapForced"],
            description="Granite area function"),

        PB.ParDouble("orgevapfrac", 0.6, 
            description="weighting for f_granitearea OrgEvapForced option" ),

        PB.ParString("f_basaltarea", "Fixed", 
            allowed_values=["Fixed", "Forced", "g3_2014_construct_from_lips"], 
            description="Basalt area function"),

        PB.ParDouble("present_basalt_area", 6.8e6, units="km^2",
            description="Present basalt area including island basalts, CFBs"),
    
        PB.ParDouble("oib_area_scaling", 2e6, units="km^2",
            description="Scaling factor to calculate ocean island basalt area from normalized DEGASS forcing"),
    )
end

function PB.register_methods!(rj::ReactionLandArea)

    vars = PB.VariableReaction[
        PB.VarPropScalar("GRAN_AREA",  "",  "Granite area normalized to present day"),
        PB.VarPropScalar("BA_AREA",  "",  "Basalt area normalized to present day"),
    ]

    if rj.pars.f_granitearea.v == "Fixed"
        # no additional Variables needed
    elseif rj.pars.f_granitearea.v == "Forced"
        push!(vars, PB.VarDepScalar("global.GRAN",     "",  "Granite area forcing normalized to present-day"))
    elseif rj.pars.f_granitearea.v == "OrgEvapForced"
        push!(vars, PB.VarDepScalar("global.GRAN",     "",  "Granite area forcing normalized to present-day"))
        push!(vars, PB.VarDepScalar("global.ORGEVAP_AREA", "",  "Contribution of silceous and shale_coal areas to overall non-basalt silicate weathering flux"))
    else
        error("$(PB.fullname(rj)) configuration error invalid f_granitearea=$(rj.pars.f_granitearea.v)")
    end

    if rj.pars.f_basaltarea.v == "Fixed"
        # no additional Variables needed
    elseif rj.pars.f_basaltarea.v == "Forced"
        push!(vars, PB.VarDepScalar("global.BA",     "",  "Basalt area forcing normalized to present-day"))
    elseif rj.pars.f_basaltarea.v == "g3_2014_construct_from_lips"
        push!(vars,
            PB.VarDepScalar("global.CFB_area", "km^2",  "Continental flood basalt area"),
            PB.VarDepScalar("global.DEGASS",   "km^2",  "Normalized degass forcing"),
            PB.VarPropScalar("oib_area",  "km^2",  "Ocean island basalt area"),
        )
    else
        error("$(PB.fullname(rj)) configuration error invalid f_basaltarea=$(rj.pars.f_basaltarea.v)")
    end  

    PB.add_method_do!(rj, do_land_area, (PB.VarList_namedtuple(vars),) )

    return nothing
end

function do_land_area(m::PB.ReactionMethod, (vars, ), cellrange::PB.AbstractCellRange, deltat)
    rj = m.reaction

    if rj.pars.f_granitearea.v == "Fixed"
        vars.GRAN_AREA[] = 1.0
    elseif rj.pars.f_granitearea.v == "Forced"
        vars.GRAN_AREA[] = vars.GRAN[]
    elseif rj.pars.f_granitearea.v == "OrgEvapForced"
        vars.GRAN_AREA[] = (1.0-rj.pars.orgevapfrac.v)*vars.GRAN[] + rj.pars.orgevapfrac.v*vars.ORGEVAP_AREA[]
    else
        error("$(PB.fullname(rj)) configuration error invalid f_granitearea=$(rj.pars.f_granitearea.v)")
    end

    if rj.pars.f_basaltarea.v == "Fixed"
        vars.BA_AREA[]  = 1.0
    elseif rj.pars.f_basaltarea.v == "Forced"
        vars.BA_AREA[]  = vars.BA[]
    elseif rj.pars.f_basaltarea.v == "g3_2014_construct_from_lips"
        vars.oib_area[] = rj.pars.oib_area_scaling.v*vars.DEGASS[]
        vars.BA_AREA[] = (vars.CFB_area[] + vars.oib_area[])/rj.pars.present_basalt_area.v
    else
        error("$(PB.fullname(rj)) configuration error invalid f_basaltarea=$(rj.pars.f_basaltarea.v)")
    end  

    return nothing
end

"Install create_reactionXXX factories when module imported"
function __init__()
    PB.add_reaction_factory(ReactionLandArea)
 
    return nothing
end

end
