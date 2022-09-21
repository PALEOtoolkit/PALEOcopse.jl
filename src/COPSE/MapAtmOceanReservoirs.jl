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
Base.@kwdef mutable struct ReactionAtmOcean_O{P} <:  PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParBool("const", false,
            description="true to provide constant value, ignoring fluxes from _sms Variable"),
    )

    norm_value::Float64 = NaN
end

function PB.register_methods!(rj::ReactionAtmOcean_O)

    # callback function to store Variable norm during setup
    function setup_callback(m, attribute_value, v, vdata)
        v.localname == "O" || error("setup_callback unexpected Variable $(PB.fullname(v))")
        if attribute_value == :norm_value
            m.reaction.norm_value = PB.value_ad(vdata[])
        end
        return nothing
    end

    if rj.pars.const[]
        O = PB.VarPropScalar(           "O",        "mol",      "atm-ocean oxygen")
        O_sms = PB.VarTargetScalar(     "O_sms",    "mol yr-1", "atm-ocean oxygen source-sinks")
        PB.add_method_setup_initialvalue_vars_default!(
            rj, [O], 
            filterfn = v->true, # force setup even though O is not a state Variable
            force_initial_norm_value=true, # force initialize as a state variable including :norm_value, even though O is not a state Variable
            setup_callback=setup_callback
        )
    else
        O = PB.VarStateExplicitScalar(  "O",        "mol",      "atm-ocean oxygen")
        O_sms = PB.VarDerivScalar(      "O_sms",    "mol yr-1", "atm-ocean oxygen source-sinks")
        PB.add_method_setup_initialvalue_vars_default!(rj, [O], setup_callback=setup_callback)  # initialise state variables 
    end
    PB.setfrozen!(rj.pars.const) # prevent modification

    # sms variable not used by us, but must appear in a method to be linked and created
    PB.add_method_do_nothing!(rj, [O_sms])

    vars = [
        O,
        PB.VarPropScalar(           "O_norm",       "",         "atm-ocean normalized"),
        PB.VarPropScalar(           "pO2atm",   "atm",      "atmospheric pO2"),
        PB.VarPropScalar(           "pO2PAL",   "",         "atmospheric pO2 normalized to present day"),
    ]

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
        PB.ParBool("const", false,
            description="true to provide constant value, ignoring fluxes from _sms Variable"),
        PB.ParType(PB.AbstractData, "CIsotope", PB.ScalarData,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable carbon isotopes and specify isotope type"),
    )

    norm_value::Float64 = NaN
end

function PB.register_methods!(rj::ReactionAtmOcean_A)

    CIsotopeType = rj.pars.CIsotope[]
  
    # callback function to store Variable norm during setup
    function setup_callback(m, attribute_value, v, vdata)
        v.localname == "A" || error("setup_callback unexpected Variable $(PB.fullname(v))")
        if attribute_value == :norm_value
            m.reaction.norm_value = PB.value_ad(PB.get_total(vdata[]))
        end
        return nothing
    end

    if rj.pars.const[]
        A = PB.VarPropScalar(           "A",          "mol",      "atm-ocean inorganic carbon (CO2 + DIC)",
            attributes=(:field_data=>CIsotopeType,))
        A_sms = PB.VarTargetScalar(     "A_sms",      "mol yr-1", "atm-ocean inorganic carbon (CO2 + DIC) source-sinks",
            attributes=(:field_data=>CIsotopeType,))
        PB.add_method_setup_initialvalue_vars_default!(
            rj, [A], 
            filterfn = v->true, # force setup even though A is not a state Variable
            force_initial_norm_value=true, # force initialize as a state variable including :norm_value, even though A is not a state Variable
            setup_callback=setup_callback
        )
    else
        A = PB.VarStateExplicitScalar(  "A",          "mol",      "atm-ocean inorganic carbon (CO2 + DIC)",
            attributes=(:field_data=>CIsotopeType,))
        A_sms = PB.VarDerivScalar(      "A_sms",      "mol yr-1", "atm-ocean inorganic carbon (CO2 + DIC) source-sinks",
            attributes=(:field_data=>CIsotopeType,))
        PB.add_method_setup_initialvalue_vars_default!(rj, [A], setup_callback=setup_callback)  # initialise state variables 
    end
    PB.setfrozen!(rj.pars.const) # prevent modification
    
    # sms variable not used by us, but must appear in a method to be linked and created
    PB.add_method_do_nothing!(rj, [A_sms])

    vars = [
        A,
        PB.VarPropScalar(         "A_norm",     "",         "atm-ocean inorganic carbon (CO2 + DIC) normalized to present day"),
        PB.VarPropScalar(         "pCO2atm","atm",      "atmospheric pCO2"),
        PB.VarPropScalar(         "pCO2PAL","",         "atmospheric pCO2 normalized to present day"),
        PB.VarPropScalar(         "phi",        "",         "atmospheric pCO2 fraction"),
    ]

    PB.add_method_do!(
        rj,
        do_AtmOcean_A, 
        (PB.VarList_namedtuple(vars),), 
    )

    if CIsotopeType <: PB.AbstractIsotopeScalar
        vars_isotope = [
            PB.VarDep(A),            
            PB.VarPropScalar("A_delta", "per mil",  "atm-ocean inorganic carbon (CO2 + DIC) delta13C"),                                    
        ]
        if rj.pars.delta_atm_ocean[]
            push!(vars_isotope,
                PB.VarDepScalar("phi", "",  "atmospheric pCO2 fraction"),
                PB.VarPropScalar("D_atmCO2_A", "per mil", "d13C fractionation between atmospheric CO2 and global DIC+CO2"),
                PB.VarPropScalar("CO2_delta", "per mil",  "atmospheric pCO2 delta 13C"),
                
                PB.VarPropScalar("D_oceanDIC_A", "per mil", "d13C fractionation between marine DIC and global DIC+CO2"),
                PB.VarPropScalar("DIC_delta", "per mil",  "marine DIC delta 13C"),
            )
        end
        if !rj.pars.fix_cisotopefrac_T[]
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
    pars,
    (vars, ),
    cellrange::PB.AbstractCellRange,
    deltat
)
    rj = m.reaction    
    
    vars.A_norm[]  = PB.get_total(vars.A[])/rj.norm_value

    if pars.f_atfrac[] == "original"
        vars.pCO2PAL[]  = vars.A_norm[]
        vars.pCO2atm[]  = vars.pCO2PAL[] * PB.Constants.k_preindCO2atm
        vars.phi[]      = 0.01614
    elseif pars.f_atfrac[] == "quadratic"
        vars.pCO2PAL[]  = vars.A_norm[]^2
        vars.pCO2atm[]  = vars.pCO2PAL[] * PB.Constants.k_preindCO2atm
        vars.phi[]      = 0.01614*vars.A_norm[]
    else
        error("unrecognized pars.f_atfrac =", pars.f_atfrac[])
    end

    return nothing
end

function do_AtmOcean_A_isotope(
    m::PB.ReactionMethod,
    pars,
    (vars, ),
    cellrange::PB.AbstractCellRange,
    deltat
)

    # fractionation of global reservoir    
    vars.A_delta[] = PB.get_delta(vars.A[])
    
    if pars.delta_atm_ocean[]
        # fractionations relative to global reservoir
        if pars.fix_cisotopefrac_T[]
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
