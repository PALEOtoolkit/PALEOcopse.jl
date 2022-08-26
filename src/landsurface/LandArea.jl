module LandArea

import PALEOboxes as PB
using PALEOboxes.DocStrings

"""
    ReactionLandArea

Calculate land areas and lithology.

Provides `GRAN_AREA` and `BA_AREA`.

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
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

    if rj.pars.f_granitearea[] == "Fixed"
        # no additional Variables needed
    elseif rj.pars.f_granitearea[] == "Forced"
        push!(vars, PB.VarDepScalar("global.GRAN",     "",  "Granite area forcing normalized to present-day"))
    elseif rj.pars.f_granitearea[] == "OrgEvapForced"
        push!(vars, PB.VarDepScalar("global.GRAN",     "",  "Granite area forcing normalized to present-day"))
        push!(vars, PB.VarDepScalar("global.ORGEVAP_AREA", "",  "Contribution of silceous and shale_coal areas to overall non-basalt silicate weathering flux"))
    else
        error("$(PB.fullname(rj)) configuration error invalid f_granitearea=$(rj.pars.f_granitearea[])")
    end

    if rj.pars.f_basaltarea[] == "Fixed"
        # no additional Variables needed
    elseif rj.pars.f_basaltarea[] == "Forced"
        push!(vars, PB.VarDepScalar("global.BA",     "",  "Basalt area forcing normalized to present-day"))
    elseif rj.pars.f_basaltarea[] == "g3_2014_construct_from_lips"
        push!(vars,
            PB.VarDepScalar("global.CFB_area", "km^2",  "Continental flood basalt area"),
            PB.VarDepScalar("global.DEGASS",   "km^2",  "Normalized degass forcing"),
            PB.VarPropScalar("oib_area",  "km^2",  "Ocean island basalt area"),
        )
    else
        error("$(PB.fullname(rj)) configuration error invalid f_basaltarea=$(rj.pars.f_basaltarea[])")
    end  

    PB.add_method_do!(rj, do_land_area, (PB.VarList_namedtuple(vars),) )

    return nothing
end

function do_land_area(m::PB.ReactionMethod, pars, (vars, ), cellrange::PB.AbstractCellRange, deltat)
    rj = m.reaction

    if pars.f_granitearea[] == "Fixed"
        vars.GRAN_AREA[] = 1.0
    elseif pars.f_granitearea[] == "Forced"
        vars.GRAN_AREA[] = vars.GRAN[]
    elseif pars.f_granitearea[] == "OrgEvapForced"
        vars.GRAN_AREA[] = (1.0-pars.orgevapfrac[])*vars.GRAN[] + pars.orgevapfrac[]*vars.ORGEVAP_AREA[]
    else
        error("$(PB.fullname(rj)) configuration error invalid f_granitearea=$(pars.f_granitearea[])")
    end

    if pars.f_basaltarea[] == "Fixed"
        vars.BA_AREA[]  = 1.0
    elseif pars.f_basaltarea[] == "Forced"
        vars.BA_AREA[]  = vars.BA[]
    elseif pars.f_basaltarea[] == "g3_2014_construct_from_lips"
        vars.oib_area[] = pars.oib_area_scaling[]*vars.DEGASS[]
        vars.BA_AREA[] = (vars.CFB_area[] + vars.oib_area[])/pars.present_basalt_area[]
    else
        error("$(PB.fullname(rj)) configuration error invalid f_basaltarea=$(pars.f_basaltarea[])")
    end  

    return nothing
end

end
