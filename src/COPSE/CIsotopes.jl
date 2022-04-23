
module CIsotopes

import PALEOboxes as PB

"""
    ReactionCIsotopes

Carbon isotope fractionation from Bergman (2004) COPSE biogeochemical model. 

Calculates global mean isotope fractionation for marine carbonate burial and land and ocean organic carbon burial
from global mean `TEMP`, `pCO2PAL` and `pO2PAL`.
"""
Base.@kwdef mutable struct ReactionCIsotopes{P} <:  PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParBool("do_D_mccb_DIC", true, 
            description="provide marine calcite burial relative to ocean DIC"),

        PB.ParBool("do_D_B_mccb_mocb", true, 
            description="provide fractionation between marine organic and calcite burial"),

        PB.ParBool("do_D_P_CO2_locb", true,
            description="provide fractionation between terrestrial organic burial and atmospheric CO2"),

        PB.ParString("f_cisotopefrac", "copse_base",
            allowed_values=["fixed", "copse_base", "copse_noO2"],
            description="fractionation calculation method"),
    )
end

function PB.register_methods!(rj::ReactionCIsotopes)

    vars = PB.VariableReaction[
        PB.VarPropScalar("land.D_eqbw_CO2", "per mil", "d13C fractionation between atmospheric CO2 and fresh water")
    ]

    if rj.pars.do_D_mccb_DIC.v
        push!(vars, PB.VarPropScalar("ocean.D_mccb_DIC", "per mil", "d13C marine calcite burial relative to ocean DIC"))
    end
    if rj.pars.do_D_B_mccb_mocb.v
        push!(vars, PB.VarPropScalar("ocean.D_B_mccb_mocb", "per mil", "d13C fractionation between marine organic and calcite burial"))
    end
    if rj.pars.do_D_P_CO2_locb.v
        push!(vars, PB.VarPropScalar("land.D_P_CO2_locb", "per mil", "d13C fractionation between terrestrial organic burial and atmospheric CO2"))
    end
             
    if rj.pars.f_cisotopefrac.v in ("copse_base", "copse_noO2")
        push!(vars,
            PB.VarDepScalar("atm.pCO2PAL", "",  "atmospheric pCO2 normalized to present day"),
            PB.VarDepScalar("TEMP", "K", "global surface temperature"),
        )
    end

    if rj.pars.f_cisotopefrac.v == "copse_base"
        push!(vars, PB.VarDepScalar("atm.pO2PAL", "",  "atmospheric pO2 normalized to present day"))
    end

    PB.add_method_do!(rj, do_CIsotopes, (PB.VarList_namedtuple(vars),))

    return nothing
end


function do_CIsotopes(
    m::PB.ReactionMethod, 
    (vars, ), 
    cellrange::PB.AbstractCellRange,
    deltat
)
    rj = m.reaction

    if      rj.pars.f_cisotopefrac.v == "fixed"
        # 15 degC, 1.0*PAL CO2 and O2
        Tkelvin = 15.0 + PB.Constants.k_CtoK
        pCO2PAL = 1.0
        pO2PAL = 1.0
    elseif  rj.pars.f_cisotopefrac.v == "copse_base"
        Tkelvin = vars.TEMP[]
        pCO2PAL = vars.pCO2PAL[]
        pO2PAL =  vars.pO2PAL[]
    elseif  rj.pars.f_cisotopefrac.v == "copse_noO2"
        Tkelvin = vars.TEMP[]
        pCO2PAL = vars.pCO2PAL[]
        pO2PAL =  1.0
    else
        error("unrecognized f_cisotopefrac=", rj.pars.f_cisotopefrac.v)
    end

    if rj.pars.do_D_mccb_DIC.v
        # marine calcite burial relative to ocean DIC
        vars.D_mccb_DIC[] = 15.10 - 4232.0/Tkelvin
    end

    if rj.pars.do_D_B_mccb_mocb.v
        # fractionation between marine organic and calcite burial
        pCO2PAL_min = 1e-3 # guard against -ve pCO2PAL when using AD
        vars.D_B_mccb_mocb[] = 33.0 -9.0/sqrt(max(pCO2PAL, pCO2PAL_min)) + 5*(pO2PAL-1)
    end

    if rj.pars.do_D_P_CO2_locb.v
        # fractionation between terrestrial organic burial and atmospheric CO2
        vars.D_P_CO2_locb[] = 19 + 5*(pO2PAL-1)
    end

    # fractionation between atmosphere and fresh (riverine) water
    # used to calculate d13C of runoff (not critical to get this right, as there is a 'short circuit' atm <-> river -> ocean <-> atm)
    vars.D_eqbw_CO2[] = 10.78-0.114*(Tkelvin - PB.Constants.k_CtoK)

    return nothing
end

# Install create_reactionXXX factories when module imported
function __init__()
    PB.add_reaction_factory(ReactionCIsotopes)
    return nothing
end

end
