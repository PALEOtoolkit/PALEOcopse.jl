# -*- coding: utf-8 -*-
module SedCrustCOPSE

import PALEOboxes as PB
using PALEOboxes.DocStrings
import PALEOcopse


"""
    ReactionSedCrustCOPSE 

COPSE Bergman(2004), COPSE Reloaded Lenton etal (2018) metamorphic and volcanic "degassing" of
sedimentary carbon and sulphur reservoirs.

Fluxes are added to flux couplers:
- `fluxSedCrusttoAOcean`:  carbon and sulphur "degassing" fluxes

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionSedCrustCOPSE{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        ##### model constants (Table 4 Bergman etal 2004, consts.hh in COPSE 5_14 C code)
        PB.ParDouble("k12_ccdeg",  6.65e12,    units="mol/yr",
            description="carbonate C degassing rate"),
        PB.ParDouble("k13_ocdeg",  1.25e12,    units="mol/yr",
            description="org C degassing rate"),

        #### COPSE S system
        PB.ParBool("enableS",          true,
            description="enable S degassing"),
        PB.ParDouble("k_pyrdeg",   0.0,        units="mol S/yr",
            description="pyrite degassing rate"),
        PB.ParDouble("k_gypdeg",   0.0,        units="mol S/yr",
            description="gypsum degassing rate"),
  
        # Organic carbon degassing
        # COPSE 5_14 C code (and Berman 2004) use "O2copsecrashprevent" which rolls
        # off organic carbon degassing at low pO2. This has a big effect at low pO2 when
        # oxidative weathering is oxygen-independent (ie Ordovician and earlier)
        PB.ParString("f_ocdeg",    "O2indep",  allowed_values=["O2indep", "O2copsecrashprevent"],
            description="roll off orgC degassing at low pO2"),

        # Isotopes
        PB.ParType(PB.AbstractData, "CIsotope", PB.ScalarData,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable carbon isotopes and specify isotope type"),
        PB.ParType(PB.AbstractData, "SIsotope", PB.ScalarData,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable sulphur isotopes and specify isotope type"),
    )

end


function PB.register_methods!(rj::ReactionSedCrustCOPSE)

    # isotope Types
    CIsotopeType = rj.pars.CIsotope.v
    if rj.pars.enableS.v
        SIsotopeType = rj.pars.SIsotope.v
    else
        SIsotopeType = PB.ScalarData
    end
        
    # state variables we use
    state_varnames = [
        ("C::$CIsotopeType",          "mol C",    "Sedimentary carbonate"),
        ("G::$CIsotopeType",          "mol C",    "Sedimentary organic carbon"),
    ]
    if rj.pars.enableS.v
        push!(state_varnames,
            ("GYP::$SIsotopeType",        "mol S",    "Sedimentary gypsum"),
            ("PYR::$SIsotopeType",        "mol S",    "Sedimentary pyrite"),
        )
    end
    vars_res, vars_sms, vars_dep_res = PB.Reservoirs.ReservoirLinksVector(
        Dict(), state_varnames
    )
    
    # dependencies
    vars_dep = PB.VarVector(
        PB.VarDepScalar,
        [
            ("global.tforce",  "yr"   , "time for external forcings"),         
            ("global.DEGASS",  "",      "degassing scaling"),
            ("global.Bforcing","",      "calcerous plankton evolution"),
            ("atm.pO2PAL",  ""   , "atmospheric oxygen normalized to present-day")
        ]
    )

    # properties we calculate
    prop_varnames = [    
        # degassing
        ("ccdeg",   "molC/yr",    "Carbonate degassing"),
        ("ocdeg",   "molC/yr",    "Organic carbon degassing"),
    ]
    if rj.pars.enableS.v
        push!(prop_varnames,
            ("gypdeg",  "molS/yr",    "Gypsum degassing"),
            ("pyrdeg",  "molS/yr",    "Pyrite degassing"),
        )
    end
    vars_prop = PB.VarVector(PB.VarPropScalar, prop_varnames)

    aocean_fluxnames = ["C::$CIsotopeType", "Redox"]
    if rj.pars.enableS.v
        push!(aocean_fluxnames, "S::$SIsotopeType")
    end
    fluxSedCrusttoAOcean = PB.Fluxes.FluxContribScalar(
        "fluxSedCrusttoAOcean.flux_", aocean_fluxnames,
        isotope_data=Dict()
    )

    PB.add_method_do!(
        rj,
        do_sed_crust_COPSE,
        (
            PB.VarList_namedtuple_fields(fluxSedCrusttoAOcean),
            PB.VarList_namedtuple([vars_res; vars_sms]),
            PB.VarList_namedtuple([vars_dep_res; vars_dep; vars_prop]),
        ),
        p=(CIsotopeType, SIsotopeType),
    )

    return nothing
end


function do_sed_crust_COPSE(
    m::PB.ReactionMethod,
    (fluxSedCrusttoAOcean, S, D),
    cellrange::PB.AbstractCellRange,
    deltat
)
    pars = m.reaction.pars
    (CIsotopeType, SIsotopeType) = m.p
    
    #%%%%%%%% calculate degassing
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    #  Inorganic carbon
    D.ccdeg[]  = pars.k12_ccdeg.v*D.DEGASS[]*D.C_norm[]*D.Bforcing[]
    ccdeg_isotope = @PB.isotope_totaldelta(CIsotopeType, D.ccdeg[], D.C_delta[])   

    # Organic carbon
    ocdeg_raw = pars.k13_ocdeg.v*D.DEGASS[]*D.G_norm[]
    if pars.f_ocdeg.v == "O2indep"
        D.ocdeg[] = ocdeg_raw
    elseif pars.f_ocdeg.v == "O2copsecrashprevent"
        # COPSE 5_14 does this (always) apparently to prevent pO2 dropping to zero ?
        # This has a big effect when pO2 dependence of oxidative weathering switched off
        D.ocdeg[] = ocdeg_raw*PALEOcopse.COPSE.copse_crash(D.pO2PAL[], "ocdeg", D.tforce[])
    else
        error("unrecogized pars.f_ocdeg ", pars.f_ocdeg.v)
    end
    
    ocdeg_isotope = @PB.isotope_totaldelta(CIsotopeType, D.ocdeg[], D.G_delta[])
  
    # update external fluxes, and state variables
    
    S.C_sms[]                   -= ccdeg_isotope    
    S.G_sms[]                   -= ocdeg_isotope  
    fluxSedCrusttoAOcean.C[]    += ccdeg_isotope + ocdeg_isotope
    # oxidant flux to AOcean
    fluxSedCrusttoAOcean.Redox[] += -D.ocdeg[]

    if pars.enableS.v
        # Sulphur
        D.pyrdeg[] = pars.k_pyrdeg.v*D.PYR_norm[]*D.DEGASS[]
        D.gypdeg[] = pars.k_gypdeg.v*D.GYP_norm[]*D.DEGASS[]
        pyrdeg_isotope = @PB.isotope_totaldelta(SIsotopeType, D.pyrdeg[], D.PYR_delta[])
        gypdeg_isotope = @PB.isotope_totaldelta(SIsotopeType, D.gypdeg[], D.GYP_delta[])
    
        # update external fluxes, and state variables
        S.GYP_sms[]                 -= gypdeg_isotope
        S.PYR_sms[]                 -= pyrdeg_isotope
        fluxSedCrusttoAOcean.S[]    += gypdeg_isotope + pyrdeg_isotope
        
        # oxidant flux to AOcean
        fluxSedCrusttoAOcean.Redox[] += - 2*D.pyrdeg[]
    end

    return nothing
end


end # module
