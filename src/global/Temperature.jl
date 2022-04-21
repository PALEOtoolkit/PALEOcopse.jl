# -*- coding: utf-8 -*-
module Temperature


import PALEOboxes as PB

"""
    ReactionGlobalTemperatureCK1992
    
Global temperature iterative calculation from SOLAR forcing, pCO2atm.
Provides a `TEMP` state variable, with either an error term for a DAE solver (if `temp_DAE = true`), 
or a restoring term for an ODE solver (if `temp_DAE = false`).

Valid range: 1e-8 < pCO2atm < 1e-2 (bar), 273.15 < TEMP < 283.15 (K)
"""
Base.@kwdef mutable struct ReactionGlobalTemperatureCK1992{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("fixed_albedo", NaN,
            description="NaN to calculate albedo, value in range (0.0 - 1.0) to fix albedo"),
        PB.ParDouble("tempcorrect", 0.194, units="K",
            description="COPSE temp correction to 15C at present day"),
        PB.ParBool("temp_DAE", false,
            description="solve TEMP as an algebraic constraint (requires DAE solver"),
    )

end

function PB.register_methods!(rj::ReactionGlobalTemperatureCK1992)
    
    vars = PB.VariableReaction[
        PB.VarDepScalar(  "atm.pCO2atm", "atm", "atmospheric CO2 partial pressure"),
        PB.VarDepScalar(  "SOLAR", "W m-2", "incident solar radiation at Earth"),
        PB.VarPropScalar(  "albedo", "", "planetary albedo"),
        PB.VarPropScalar(  "Teff", "K", "black body effective temperature"),
        PB.VarPropScalar(  "Tgreenhouse", "K", ""),
    ]
       
    if rj.pars.temp_DAE.v
        push!(vars, PB.VarStateScalar("TEMP", "K", "global surface temperature",
                        attributes=(:initial_value=>PB.Constants.k_CtoK+15.0, ))
        )
        push!(vars, PB.VarConstraintScalar("TEMP_constraint", "K", "algebraic constraint (==0)",
                        attributes=(:norm_value=>1.0, ))
        )
    else
        push!(vars, PB.VarStateExplicitScalar("TEMP", "K", "global surface temperature",
                        attributes=(:initial_value=>PB.Constants.k_CtoK+15.0, ))
        )
        push!(vars, PB.VarDerivScalar("TEMP_sms", "K yr-1", "rate of change of global surface temperature"))
    end
       
    PB.add_method_setup_initialvalue_vars_default!(rj, vars)  # initialise state variables 

    PB.add_method_do!(rj, do_GlobalTemperatureCK1992, (PB.VarList_namedtuple(vars),))

    PB.add_method_initialize_zero_vars_default!(rj)

    return nothing
end


# Calculate iterative improvement to temperature estimate
function do_GlobalTemperatureCK1992(m::PB.ReactionMethod, (vars, ), cellrange::PB.AbstractCellRange, deltat)

    rj = m.reaction

    newT, vars.albedo[], vars.Tgreenhouse[], vars.Teff[] = copse_TempCK1992(
        vars.SOLAR[], vars.pCO2atm[], vars.TEMP[], fixed_albedo=rj.pars.fixed_albedo.v
    )
    newT += rj.pars.tempcorrect.v

    # algebraic constraint or iterative improvement to temperature estimate
    if rj.pars.temp_DAE.v
        vars.TEMP_constraint[] = newT - vars.TEMP[]
    else
        vars.TEMP_sms[] = newT - vars.TEMP[]
    end
    
    return nothing
end


"""
    copse_TempCK1992( luminosity, pCO2atm, oldT [, albedo] ) -> (newT, albedo, Tgh, Teff)

Caldeira & Kasting (1992) global temperature function.

# Arguments
- `luminosity`: solar luminosity, W/m^2 (currently 1368.0)
- `pCO2atm`: (bar) atmospheric pCO2 (pre-industrial 280e-6, valid range 1e-8 - 1e-2)
- `oldT`:   (K) estimated temperature (valid range 273.15 - 283.15)
- `albedo`: (optional) fix planetary albedo, default is variable albedo

# Returns
- `newT`: (K) improved temperature estimate
- `albedo`: planetary albedo
- `Tgh`: (K) greenhouse contribution
- `Teff`: (K) black-body effective temperature

Caldeira, K., & Kasting, J. F. (1992). The life span of the biosphere revisited. Nature, 360(6406), 721–3. doi:10.1038/360721a0
"""
function copse_TempCK1992( luminosity, pCO2atm, oldT; fixed_albedo=NaN)

    CK_sigma = 5.67e-8  # Stefan-Boltzmann
    CK_a0    = 1.4891
    CK_a1    = -0.0065979
    CK_a2    = 8.567e-6
    CK_1 = 815.17
    CK_2 = 4.895e+7
    CK_3 = -3.9787e+5
    CK_4 = -6.7084
    CK_5 = 73.221
    CK_6 = -30882.0

    if isnan(fixed_albedo)        
        albedo = CK_a0 + CK_a1 * oldT + CK_a2 * oldT^2
    else
        albedo = one(oldT)*fixed_albedo  # convert to type of oldT for type stability with AD
    end

    # calculate the effective black body temperature
    Teff =((1.0-albedo) * luminosity / (4.0 * CK_sigma))^0.25

    # calculate the greenhouse effect */
    pCO2atm_min = 300e-9  # guard against -ve pCO2 when using AD
    psi = log10(max(pCO2atm, pCO2atm_min))
    Tgh = (CK_1 + CK_2 / oldT^2 + CK_3 / oldT + CK_4 / psi^2
            + CK_5 / psi + CK_6 / (psi * oldT))

    # Final T is black body T + greenhouse effect

    newT = Teff + Tgh

    # BM Matlab code - identical to above
    # tempcorrect = 0.194; %%%% COPSE temp correction to 15C at present day
    # TEMP = (       (  ((1- (1.4891 - 0.0065979*S.temp + (8.567e-6)*(S.temp^2)  )  )*SOLAR)/(4*5.67e-8)  )^0.25...
    #     +  815.17  + (4.895e7)*(S.temp^-2) -  (3.9787e5)*(S.temp^-1)...
    #     -6.7084*((log10(  CO2atm  ))^-2) + 73.221*((log10(  CO2atm  ))^-1) -30882*(S.temp^-1)*((log10(  CO2atm ))^-1)     ) + tempcorrect ;

    # CK1992 fixed albedo
    # TEMP = (       (  ((1- ( ALBEDO )  )*D.SOLAR)/(4*5.67e-8)  )^0.25     +  815.17  + (4.895e7)*(S.temp^-2) -  (3.9787e5)*(S.temp^-1)  -6.7084*((log10(  D.pCO2atm  ))^-2) + 73.221*((log10(  D.pCO2atm  ))^-1) -30882*(S.temp^-1)*((log10(  D.pCO2atm ))^-1)     ) ;

    return (newT, albedo, Tgh, Teff)
end



"""
    ReactionGlobalTemperatureBerner
    
Global temperature calculation from pCO2, implicit solar luminosity (from forcing time tforce).
"""
Base.@kwdef mutable struct ReactionGlobalTemperatureBerner{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("k_c",   4.328,  description="climate sensitivity to CO2"),
        PB.ParDouble("k_l",   7.4,    description="rate of change of solar luminosity"),
    )

end

function PB.register_methods!(rj::ReactionGlobalTemperatureBerner)
    vars = PB.VariableReaction[
        PB.VarPropScalar( "TEMP", "K", "global surface temperature")
        PB.VarDepScalar(  "atm.pCO2PAL", "", "atmospheric CO2 partial pressure relative to pre-industrial")
        PB.VarDepScalar(  "tforce", "yr",  "historical time at which to apply forcings, present = 0 yr")
    ]
    PB.add_method_do!(rj, do_GlobalTemperatureBerner, (PB.VarList_namedtuple(vars),))

    return nothing
end

function do_GlobalTemperatureBerner(m::PB.ReactionMethod, (vars, ), cellrange::PB.AbstractCellRange, deltat)
   
    rj = m.reaction

    pCO2PAL = clamp(vars.pCO2PAL[], 1e-3, 1e3) # guess at pCO2PAL limit
    pCO2PAL == vars.pCO2PAL[] ||
        @warn "ReactionGlobalTemperatureBerner pCO2PAL $(vars.pCO2PAL[]) out of range, limiting to $pCO2PAL"

    # TL note that this uses its own luminosity (implicitly)
    TEMP = (PB.Constants.k_CtoK + 15.0 +
            rj.pars.k_c.v*log(pCO2PAL) +
            rj.pars.k_l.v*vars.tforce[]/570e6)
    
    vars.TEMP[] = clamp(TEMP, PB.Constants.k_CtoK - 50.0, PB.Constants.k_CtoK + 50.0 )  # guess at TEMP limit
    vars.TEMP[] == TEMP ||
        @warn "ReactionGlobalTemperatureBerner TEMP $TEMP out of range, limiting to $(vars.TEMP[])"

    return nothing
end


"Install create_reactionXXX factories when module imported"
function __init__()
    PB.add_reaction_factory(ReactionGlobalTemperatureCK1992)
    PB.add_reaction_factory(ReactionGlobalTemperatureBerner)
    return nothing
end

end # module
