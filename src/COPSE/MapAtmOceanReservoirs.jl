# -*- coding: utf-8 -*-
module MapAtmOceanReservoirs

import PALEOboxes as PB
using PALEOboxes.DocStrings

"""
    ReactionAtmOcean_O

Combined atmosphere + ocean scalar O (oxygen) reservoir for use with COPSE model.

Provides state variable `O` (O2 in moles), state variable  time derivative `O_sms` (mol yr-1) and 
`O_norm` (`O` normalized by `O:norm_value` attribute, which should be set in .yaml config
file to present-day atmospheric value).

Also calculates additional quantities `pO2atm` (partial pressure in bar), and `pO2PAL` (equal to `O_norm`), which should usually
be relinked to the `atm` Domain in the config file. Assumes ocean oxygen is a neglible fraction of the total.

# Parameters
$(PARS)

# Methods and Variables
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionAtmOcean_O <:  PB.AbstractReaction
    base::PB.ReactionBase

    norm_value::Float64 = NaN
end

function PB.register_methods!(rj::ReactionAtmOcean_O)

    vars = [
        PB.VarStateExplicitScalar(  "O",            "mol",      "atm-ocean oxygen"),
        PB.VarDerivScalar(          "O_sms",        "mol yr-1", "atm-ocean oxygen source-sinks"),
        PB.VarPropScalar(           "O_norm",       "",         "atm-ocean normalized"),
        PB.VarPropScalar(           "pO2atm",   "atm",      "atmospheric pO2"),
        PB.VarPropScalar(           "pO2PAL",   "",         "atmospheric pO2 normalized to present day"),
    ]

    # callback function to store Variable norm during setup
    function setup_callback(m, attribute_value, v, vdata)
        v.localname == "O" || error("setup_callback unexpected Variable $(PB.fullname(v))")
        if attribute_value == :norm_value
            m.reaction.norm_value = PB.value_ad(vdata[])
        end
        return nothing
    end

    PB.add_method_setup_initialvalue_vars_default!(rj, vars, setup_callback=setup_callback)  # initialise state variables 

    PB.add_method_do!(
        rj,
        do_AtmOcean_O, 
        (PB.VarList_namedtuple(vars),), 
    )

    PB.add_method_initialize_zero_vars_default!(rj)

    return nothing
end


function do_AtmOcean_O(
    m::PB.ReactionMethod,
    (vars, ), 
    cellrange::PB.AbstractCellRange,
    deltat
)
    rj = m.reaction

    vars.O_norm[]  = vars.O[]/rj.norm_value

    vars.pO2PAL[] = vars.O_norm[]
    
    vars.pO2atm[] = vars.O[] / PB.Constants.k_moles1atm

    return nothing
end



"""
    ReactionAtmOcean_A

Atmosphere-ocean inorganic carbon `A` reservoir, for use with the COPSE model.

Provides state variable `A` (total C in atmospheric CO2 and ocean DIC, in moles), 
state variable  time derivative `A_sms` (mol yr-1) and  `A_norm` (`A` normalized by `A:norm_value` attribute,
which should be set in .yaml config file to the pre-industrial value of 3.193e18 mol).

Also calculates partitioning into atmosphere and ocean (atm-ocean fraction `phi`) and 
atmospheric `pCO2atm` (partial pressure in bar), and `pCO2PAL`, which should usually
be relinked to the `atm` Domain in the config file.

The calculation of `phi` atm-ocean fraction is set by Parameter `f_atfrac`:
- `original`:  fixed `phi`, as used in the original COPSE model [Bergman2004](@cite).
- `quadratic`: `phi` proportional to `A`, approximating the behaviour of the carbonate system assuming carbonate saturation
  hence `[CO3--]` carbonate ion concentration is approximately constant.

If Parameter `delta_atm_ocean == true`, also calculates atmosphere CO2 and ocean DIC isotopic fractionation.

# Parameters
$(PARS)

# Methods and Variables
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionAtmOcean_A{P} <:  PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParString("f_atfrac", "original",
            allowed_values=["original", "quadratic"],
            description="atm-ocean partitioning"),
        PB.ParBool("delta_atm_ocean", true, 
            description="calculate d13CO2, d13DIC relative to A"),
        PB.ParBool("fix_cisotopefrac_T", false,
            description="remove temperature dependence of d13CO2, d13DIC relative to A"),
        PB.ParType(PB.AbstractData, "CIsotope", PB.ScalarData,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable carbon isotopes and specify isotope type"),
    )

    norm_value::Float64 = NaN
end

function PB.register_methods!(rj::ReactionAtmOcean_A)

    CIsotopeType = rj.pars.CIsotope.v
  
    vars = [
        PB.VarStateExplicitScalar("A",          "mol",      "atm-ocean inorganic carbon (CO2 + DIC)",
            attributes=(:field_data=>CIsotopeType,)),
        PB.VarDerivScalar(        "A_sms",      "mol yr-1", "atm-ocean inorganic carbon (CO2 + DIC) source-sinks",
            attributes=(:field_data=>CIsotopeType,)),
        PB.VarPropScalar(         "A_norm",     "",         "atm-ocean inorganic carbon (CO2 + DIC) normalized to present day"),
        PB.VarPropScalar(         "pCO2atm","atm",      "atmospheric pCO2"),
        PB.VarPropScalar(         "pCO2PAL","",         "atmospheric pCO2 normalized to present day"),
        PB.VarPropScalar(         "phi",        "",         "atmospheric pCO2 fraction"),
    ]

    # callback function to store Variable norm during setup
    function setup_callback(m, attribute_value, v, vdata)
        v.localname == "A" || error("setup_callback unexpected Variable $(PB.fullname(v))")
        if attribute_value == :norm_value
            m.reaction.norm_value = PB.value_ad(PB.get_total(vdata[]))
        end
        return nothing
    end

    PB.add_method_setup_initialvalue_vars_default!(rj, vars, setup_callback=setup_callback)  # initialise state variables 

    norm_value = NaN # will be updated in prepare_
    PB.add_method_do!(
        rj,
        do_AtmOcean_A, 
        (PB.VarList_namedtuple(vars),), 
    )


    if CIsotopeType <: PB.AbstractIsotopeScalar
        vars_isotope = [
            PB.VarDepScalar("A",          "mol",      "atm-ocean inorganic carbon (CO2 + DIC)",
                attributes=(:field_data=>CIsotopeType,)),
            
            PB.VarPropScalar("A_delta", "per mil",  "atm-ocean inorganic carbon (CO2 + DIC) delta13C"),                                    
        ]
        if rj.pars.delta_atm_ocean.v
            push!(vars_isotope,
                PB.VarDepScalar("phi", "",  "atmospheric pCO2 fraction"),
                PB.VarPropScalar("D_atmCO2_A", "per mil", "d13C fractionation between atmospheric CO2 and global DIC+CO2"),
                PB.VarPropScalar("CO2_delta", "per mil",  "atmospheric pCO2 delta 13C"),
                
                PB.VarPropScalar("D_oceanDIC_A", "per mil", "d13C fractionation between marine DIC and global DIC+CO2"),
                PB.VarPropScalar("DIC_delta", "per mil",  "marine DIC delta 13C"),
            )
        end
        if !rj.pars.fix_cisotopefrac_T.v
            push!(vars_isotope, PB.VarDepScalar("TEMP", "K", "global surface temperature"))
        end
        PB.add_method_do!(
            rj,
            do_AtmOcean_A_isotope, 
            (PB.VarList_namedtuple(vars_isotope),), 
        )
    end

    PB.add_method_initialize_zero_vars_default!(rj)
    
    return nothing
end

function do_AtmOcean_A(
    m::PB.ReactionMethod,
    (vars, ),
    cellrange::PB.AbstractCellRange,
    deltat
)
    rj = m.reaction    
    
    vars.A_norm[]  = PB.get_total(vars.A[])/rj.norm_value

    if rj.pars.f_atfrac.v == "original"
        vars.pCO2PAL[]  = vars.A_norm[]
        vars.pCO2atm[]  = vars.pCO2PAL[] * PB.Constants.k_preindCO2atm
        vars.phi[]      = 0.01614
    elseif rj.pars.f_atfrac.v == "quadratic"
        vars.pCO2PAL[]  = vars.A_norm[]^2
        vars.pCO2atm[]  = vars.pCO2PAL[] * PB.Constants.k_preindCO2atm
        vars.phi[]      = 0.01614*vars.A_norm[]
    else
        error("unrecognized pars.f_atfrac.v=", rj.pars.f_atfrac.v)
    end

    return nothing
end

function do_AtmOcean_A_isotope(
    m::PB.ReactionMethod,
    (vars, ),
    cellrange::PB.AbstractCellRange,
    deltat
)
    rj = m.reaction

    # fractionation of global reservoir    
    vars.A_delta[] = PB.get_delta(vars.A[])
    
    if rj.pars.delta_atm_ocean.v
        # fractionations relative to global reservoir
        if rj.pars.fix_cisotopefrac_T.v
            Tkelvin = 15.0 + PB.Constants.k_CtoK
        else
            Tkelvin = vars.TEMP[]
        end

        # ocean total dissolved CO2 (DIC) relative to atm-ocean A
        vars.D_oceanDIC_A[] = vars.phi[]*(9483.0/Tkelvin-23.89)
        vars.DIC_delta[] = vars.A_delta[] + vars.D_oceanDIC_A[]
        
        # atmosphere CO2 relative to atm-ocean A
        vars.D_atmCO2_A[] = (vars.phi[]-1.0)*(9483.0/Tkelvin-23.89)
        vars.CO2_delta[] = vars.A_delta[] + vars.D_atmCO2_A[]
    end

end

end
