module CarbBurial

import PALEOboxes as PB
using PALEOboxes.DocStrings

"""
    ReactionCarbBurialAlk

Marine carbonate burial from instantaneous alkalinity balance

Generates burial flux `fluxOceanBurial.flux_Ccarb = 0.5 * flux_TAlk`, adds reservoirs fluxes to `DIC_sms` and
optionally `CAL_sms`.

Fluxes are added to flux couplers:
- `fluxOceanBurial`: ocean burial fluxes

Default linking for `flux_TAlk` is `flux_TAlk -> fluxRtoOcean.flux_TAlk`, ie carbonate burial sink for riverine alkalinity input.

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionCarbBurialAlk{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        # Isotopes
        PB.ParType(PB.AbstractData, "CIsotope", PB.ScalarData,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable carbon isotopes and specify isotope type"),        
    )

end


function PB.register_methods!(rj::ReactionCarbBurialAlk)

    # isotope Types
    CIsotopeType = rj.pars.CIsotope[]

    # dependencies required in do_react
    vars = PB.VariableReaction[      
        # Special-case a dependency on (riverine) alk flux for carbonate burial
        PB.VarDepScalar("flux_TAlk"=>"fluxRtoOcean.flux_TAlk",  "mol yr-1"   , "riverine alkalinity flux"),  
        # burial
        PB.VarPropScalar("mccb",   "molC/yr",     "Carbonate burial"),
        PB.VarContribScalar("ocean.oceanfloor.DIC_sms", "molC/yr", "DIC reservoir source - sink",
            attributes=(:field_data=>CIsotopeType,)),
        PB.VarContribScalar("(ocean.oceanfloor.CAL_sms)", "mol Ca/yr", "Calcium reservoir source - sink"),    
        PB.VarContribScalar("fluxOceanBurial.flux_Ccarb", "molC/yr", "carbonate burial",
            attributes=(:field_data=>CIsotopeType,)),
    ]

    if CIsotopeType <: PB.AbstractIsotopeScalar
        append!(vars, [
                PB.VarDepScalar("ocean.oceanfloor.DIC_delta",  "per mil",  "d13C ocean DIC"),
                PB.VarDepScalar("ocean.D_mccb_DIC",  "per mil",  "D13C marine calcite burial relative to ocean DIC"),
                # only for diagnostic output
                PB.VarDepScalar("ocean.D_oceanDIC_A", "per mil", "d13C fractionation between marine DIC and global DIC+CO2"),
                # C isotopes
                PB.VarPropScalar("ocean.mccb_delta", "per mil", "D13C fractionation of marine calcium carbonate burial"),   
                PB.VarPropScalar("ocean.D_mccb_A",    "per mil", "d13C fractionation between marine carbonate burial and global DIC+CO2"),    
            ]
        )
    end

    PB.add_method_do!(
        rj,
        do_carb_burial,
        (
            PB.VarList_namedtuple(vars),
        ),
        p=CIsotopeType
    )

    return rj
end


# Calculate rates
function do_carb_burial(m::PB.ReactionMethod,  (vars, ), cellrange::PB.AbstractCellRange, deltat)
   
    CIsotopeType = m.p

    if CIsotopeType <: PB.AbstractIsotopeScalar
        # delta of marine carbonate burial
        vars.mccb_delta[] = vars.DIC_delta[] + vars.D_mccb_DIC[]
       
        # for diagnostic output only
        vars.D_mccb_A[] = vars.D_oceanDIC_A[] + vars.D_mccb_DIC[]
    end
 
    vars.mccb[]                = 0.5*vars.flux_TAlk[]  # alkalinity balance
    mccb_isotope            = @PB.isotope_totaldelta(CIsotopeType, vars.mccb[], vars.mccb_delta[])
    
    PB.add_if_available(vars.CAL_sms, -vars.mccb[])
    vars.DIC_sms[]             += -mccb_isotope
    vars.flux_Ccarb[] += mccb_isotope

    return nothing
end

end