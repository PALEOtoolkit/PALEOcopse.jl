# -*- coding: utf-8 -*-
module OceanCOPSE


import PALEOboxes as PB
using PALEOboxes.DocStrings

import PALEOcopse

"""
    ReactionMarineBiotaCOPSE

COPSE Bergman(2004), COPSE Reloaded Lenton etal (2018) 0D marine biota

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionMarineBiotaCOPSE{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        ##### model constants (Table 4 Bergman etal 2004, consts.hh in COPSE 5_14 C code)
        PB.ParDouble("k1_oxfrac",  0.86,       units="",       description="initial oxic fraction"),
        PB.ParDouble("k3_nfix",    8.75e12,    units="mol/yr", description="nitrogen fixation"),
        PB.ParDouble("k4_denit",   4.3e12,     units="mol/yr", description="denitrification"),

        ##########################################################################
        ####### Options controlling functional forms
        ##########################################################################

        # Marine N cycle
        PB.ParBool("f_ncycle", true,
            description="enable nitrogen cycle"),
        PB.ParDouble("f_nfix_power", 2.0,      units="",       
            description="nitrogen fixation power-law dependence on nitrogen deficit wrt Redfield"),
        PB.ParString("f_nfix_nreplete","Off",  allowed_values=["Off", "Sign"],
            description="nitrogen fixation when nitrogen excess wrt Redfield (only for bug compatitibility with COPSE 5_14 C code, which has Sign)"),
        PB.ParString("f_denit",     "original", allowed_values=["original", "new"],
            description="functional form for denitrification"),

        # Marine ecosystem
        PB.ParDouble("newp0",      225.956),

        PB.ParString("f_anoxia",    "original", allowed_values=["original", "newanoxia"],
            description="functional form for marine anoxia function"),
        PB.ParDouble("k_logistic",  12.0,
            description="slope of logistic anoxia function"),
        PB.ParDouble("k_uptake",    0.5,
            description="efficiency of nutrient uptake in anoxia function"),      
    )

end


function PB.register_methods!(rj::ReactionMarineBiotaCOPSE)

    state_varnames = [
        ("P",                   "mol",    "Marine phosphorus"),         
    ]   
    if rj.pars.f_ncycle[]
        push!(state_varnames, ("N",               "mol",    "Marine nitrogen"))
    end
    vars_res, vars_sms, vars_dep_res = PB.Reservoirs.ReservoirLinksVector(
        Dict(), state_varnames
    )

    # dependencies required in do_react
    vars_dep_react = PB.VarVector(
        PB.VarDepScalar,
        [ 
            ("global.tforce",    "yr",      "epoch for model forcing"),       
            ("atm.pO2PAL",    "",      "atmospheric pO2 normalized to present day"),
        ]
    )
            
    # Properties we calculate in do_react
    varnames_prop_react = [
        # copse_marinebiota
        ("Pconc",   "molP/kgsw",  "marine P conc"),
        ("newp",    "",           "marine new production"),
        ("newp_relative",  "",    "marine new production relative to present day"),
        ("ANOX",    "",           "marine anoxic fraction"),
        ("OX_relative",    "",        "marine oxic fraction relative to present day fraction 'k1_oxfrac'"),
    ]

    if rj.pars.f_ncycle[]
        push!(varnames_prop_react,
            ("Nconc",   "molN/kgsw",  "marine N conc"),
            ("nfix",    "molN/yr",    "marine nitrogen fixation"),
            ("denit",   "molN/yr",    "marine denitrification"),
        )
    end
    vars_prop_react = PB.VarVector(PB.VarPropScalar, varnames_prop_react)

    # Now assemble these lists into properties and dependencies for each method,
    # following the Matlab COPSE convention that:
    #   S - includes reservoirs (a VarDep) and reservoir_sms (a VarContrib)
    #   D - includes dependencies (VarDep) and properties to calculate (VarProp)

    S_react = [vars_res; vars_sms]
    # Properties we calculated in do_isotopes are dependencies for do_react
    D_react = [vars_dep_res; vars_dep_react; vars_prop_react]
    PB.add_method_do!(
        rj,
        do_marinebiota,
        (
            PB.VarList_namedtuple(S_react), 
            PB.VarList_namedtuple(D_react),
        ),
    )

    return rj
end


# Calculate rates
function do_marinebiota(m::PB.ReactionMethod, pars, (S, D), cellrange::PB.AbstractCellRange, deltat)
   
    # convert marine nutrient reservoir moles to micromoles/kg concentration
    D.Pconc[] = D.P_norm[] * 2.2  
    #  clunky way to get normalization values
    # fails as Reaction doesn't see norm_value for variable it is linked to
    # P0 = S.P.norm_value
    # N0 = S.N.norm_value
    P0 = S.P[]/D.P_norm[]

    if pars.f_ncycle[]
        D.Nconc[] = D.N_norm[] * 30.9 ;
        N0 = S.N[]/D.N_norm[]
        D.newp[] = 117.0 * min(D.Nconc[]/16.0, D.Pconc[])
    else
        D.newp[] = 117.0 * D.Pconc[]
    end
    D.newp_relative[] = D.newp[]/pars.newp0[]

    #%%%%%% OCEAN ANOXIC FRACTION
    if      pars.f_anoxia[] == "original"
        D.ANOX[] = max( 1.0 - pars.k1_oxfrac[]*D.pO2PAL[]/D.newp_relative[], 0.0 )
    elseif  pars.f_anoxia[] == "newanoxia"
        D.ANOX[] = 1/(1 + exp(-pars.k_logistic[]*(pars.k_uptake[]*D.newp_relative[]-D.pO2PAL[]) ) )
    else
        error("unknown f_anoxia ", pars.f_anoxia[])
    end
    D.OX_relative[] = (1.0 - D.ANOX[])/pars.k1_oxfrac[]   

    #%%%%% nitrogen cycle
    if pars.f_ncycle[]
        if (S.N[]/16.0) < S.P[] 
            D.nfix[] = pars.k3_nfix[] *( ( S.P[] - S.N[]/16.0 ) / ( P0  - N0/16.0) )^pars.f_nfix_power[]
        else
            if   pars.f_nfix_nreplete[] == "Off" # Surely more defensible ?
                D.nfix[] = 0.0
            elseif pars.f_nfix_nreplete[] == "Sign" # SD - COPSE 5_14 C code has this (?!)
                D.nfix[] = pars.k3_nfix[] *( - ( S.P[] - S.N[]/16.0 ) / ( P0  - N0/16.0)  )^pars.f_nfix_power[]
                print("COPSE_equations -ve nfix (check pars.f_nfix_nreplete) tmodel ", D.tforce[])
            else
                error("unrecognized f_nfix_nreplete ", pars.f_nfix_nreplete[])
            end
        end
        # Denitrification 
        if pars.f_denit[] == "original" # NB: COPSE 5_14 uses copse_crash to limit at low N
            D.denit[] = (pars.k4_denit[] * (1.0 + D.ANOX[] / (1.0 - pars.k1_oxfrac[]) ) *
                        PALEOcopse.COPSE.copse_crash(D.N_norm[], "denit", D.tforce[]))
        elseif  pars.f_denit[] == "new" # introduce dependency on [NO3] throughout
            D.denit[] = pars.k4_denit[] * (1.0 + D.ANOX[] / (1.0 - pars.k1_oxfrac[])) * D.N_norm[]
        else
            error("unknown f_denit ", pars.f_denit[])
        end

        S.N_sms[] += D.nfix[] - D.denit[]
    end

    return nothing
end


"""
    ReactionOceanBurialCOPSE

COPSE Bergman(2004), COPSE Reloaded Lenton etal (2018) 0D ocean N, P, Corg, S burial 

Fluxes are added to flux couplers:
- `fluxOceanBurial`: ocean burial fluxes

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionOceanBurialCOPSE{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        ##### model constants (Table 4 Bergman etal 2004, consts.hh in COPSE 5_14 C code)
        PB.ParDouble("k2_mocb",    4.5e12,     units="mol/yr", description="ocean organic carbon burial"),
        PB.ParDouble("k6_fepb",    6e9,        units="mol/yr", description="Fe-P burial"),
        PB.ParDouble("k7_capb",    1.5e10,     units="mol/yr", description="Ca-P burial"),
               
        PB.ParBool("enableS",          true,    description="enable S burial"),
        PB.ParDouble("k_mpsb",     0.53e12,    units="mol S/yr",description="pyrite burial"),
        PB.ParDouble("k_mgsb",     1e12,       units="mol S/yr",description="gypsum burial"),


        ##########################################################################
        ####### Options controlling functional forms
        ##########################################################################

        # Marine N cycle
        PB.ParBool("f_ncycle", true,
            description="enable nitrogen cycle"),
    
        # Marine CPN ratio
        PB.ParString("f_CPsea",    "Fixed",    allowed_values=["Fixed", "VCI"],
            description="Functional form of CPsea ratio"),
        PB.ParDouble("CPsea0",      250.0,
            description="for f_CPsea:\"Fixed\""),
        # Van Cappellen & Ingall (f_CPsea:"VCI") marine C/P burial ratio consts
        # These are Redfield Revisited 1 steady-state values
        PB.ParDouble("f_CPsea_VCI_oxic",217.0,
            description="Van Cappellen & Ingall (f_CPsea:\"VCI\") marine C/P burial ratio consts"),
        PB.ParDouble("f_CPsea_VCI_anoxic",4340.0,
            description="Van Cappellen & Ingall (f_CPsea:\"VCI\") marine C/P burial ratio consts"),
        #f_CPsea_VCI_oxic : 200.0     # Original VC&I values
        #f_CPsea_VCI_anoxic : 4000.0
        PB.ParString("f_capb",  "original",     allowed_values=["original", "redox"],
            description="functional form of marine carbonate-associated P burial"),
        PB.ParString("f_fepb",  "original",     allowed_values=["original", "Dforced", "sfw", "pdep"],
            description="functional form of marine iron-associated P burial"),

        PB.ParDouble("CNsea0",      37.5,
            description="Always fixed"),

        PB.ParString("f_mocb",      "original", allowed_values=["original", "Uforced", "O2dep", "both"],
            description="functinal form of marine organic carbon burial"),
        PB.ParDouble("f_mocb_b",    2.0,        
            description="marine organic carbon burial power-law dependency on new production"),

        #Marine pyrite sulphur burial dependency on oxygen
        PB.ParString("f_pyrburial", "copse_O2", allowed_values=["copse_O2", "copse_noO2"],
            description="Marine pyrite sulphur burial dependency on oxygen"),
        PB.ParBool("SRedoxAlk",     false,
            description="true to include +TAlk from pyrite burial"),
        # S isotope fractionation calculation
        PB.ParString("f_sisotopefrac","fixed",  allowed_values=["fixed", "copse_O2"],
            description="S isotope fractionation calculation."),

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

function PB.set_model_geometry(rj::ReactionOceanBurialCOPSE, model::PB.Model)
    # TODO - not the right place for this ? 
    # define a 1 cell grid
    rj.domain.grid = PB.Grids.UnstructuredVectorGrid(ncells=1)
    PB.Grids.set_subdomain!(rj.domain.grid, "oceansurface", PB.Grids.BoundarySubdomain([1]), true)
    PB.Grids.set_subdomain!(rj.domain.grid, "oceanfloor",   PB.Grids.BoundarySubdomain([1]), true)

    return nothing
end

function PB.register_methods!(rj::ReactionOceanBurialCOPSE)

    # isotope Types
    CIsotopeType = rj.pars.CIsotope[]
    if rj.pars.enableS[]
        SIsotopeType = rj.pars.SIsotope[]
    else
        SIsotopeType = PB.ScalarData
    end

    state_varnames = [          
        ("O",                   "mol",   "Atm-ocean oxygen (mol O2)"),
        ("P",                   "mol",    "Marine phosphorus"),
        ("(DIC::$CIsotopeType)",    "mol",    "ocean inorganic carbon"),    
        ("(CAL)",               "mol",   "Marine calcium"),
    ]
    if rj.pars.enableS[]
        push!(state_varnames, ("S::$SIsotopeType",      "mol",    "Marine sulphate"))
    end
    vars_res, vars_sms, vars_dep_res = PB.Reservoirs.ReservoirLinksVector(
        Dict(), state_varnames
    )

    # additional fluxes where we don't require the values of the state variable 
    if rj.pars.f_ncycle[]
        push!(vars_sms, PB.VarContribScalar("N_sms", "mol yr-1", "Marine nitrogen source - sink"))
    end
    if rj.pars.enableS[] && rj.pars.SRedoxAlk[]
        push!(vars_sms, PB.VarContribScalar("TAlk_sms", "mol yr-1", "Marine total alkalinity source - sink"))
    end

    # dependencies required in do_react
    vars_dep_react = PB.VarVector(
        PB.VarDepScalar,
        [ 
            ("(D_mccb_DIC)",  "per mil",  "D13C marine calcite burial relative to ocean DIC"),
            ("(D_B_mccb_mocb)","per mil", "D13C fractionation between marine organic and calcite burial"),
            # only for diagnostic output
            ("(D_oceanDIC_A)", "per mil", "d13C fractionation between marine DIC and global DIC+CO2"),

            ("global.tforce",    "yr",      "epoch for model forcing"),
            ("global.TEMP",    "K",      "global temperature"),
            ("(global.UPLIFT)",    "",      "UPLIFT forcing"),
            ("(global.DEGASS)",    "",      "DEGASS forcing"),            
            ("atm.pO2PAL",    "",      "atmospheric pO2 normalized to present day"),
            ("(oceanfloor.sfw_relative)",    "",      "seafloor weathering relative to present"),
            # from copse_marinebiota
            ("newp_relative",  "",        "normalized marine new production"),
            ("ANOX",    "",           "marine anoxic fraction"),
            ("OX_relative",  "",           "marine oxic fraction relative to present day"),
        ]
    )
            
    # Properties we calculate in do_react
    varnames_prop_react = [
        # C isotopes
        ("mocb_delta", "per mil", "D13C fractionation of marine organic carbon burial"),

        # burial
        ("CPsea",   "",           "marine C:P ratio"),
        ("mocb",   "molC/yr",     "Marine organic carbon burial"),       
    ]
    if rj.pars.enableS[]
        push!(varnames_prop_react, ("D_mpsb", "per mil",     "D34S fractionation pyrite - marine sulphate"))
    end
    if rj.pars.f_ncycle[]
        push!(varnames_prop_react,
            ("monb",   "molN/yr",     "Marine organic N burial"),
        )
    end
    vars_prop_react = PB.VarVector(PB.VarPropScalar, varnames_prop_react)

    # Define flux coupler (a NamedTuple of Variables)
    burial_fluxnames = [   
        "Corg::$CIsotopeType", 
        "Porg", "Pauth", "PFe", "P", 
    ]
    if rj.pars.enableS[]
        push!(burial_fluxnames, "GYP::$SIsotopeType", "PYR::$SIsotopeType")
    end

    fluxOceanBurial = PB.Fluxes.FluxContrib( # Scalar ?
        "fluxOceanBurial.flux_",
        burial_fluxnames,
        isotope_data=Dict(),
    ) 

    # Now assemble these lists into properties and dependencies for each method,
    # following the Matlab COPSE convention that:
    #   S - includes reservoirs (a VarDep) and reservoir_sms (a VarContrib)
    #   D - includes dependencies (VarDep) and properties to calculate (VarProp)

    S_react = [vars_res; vars_sms]
    # Properties we calculated in do_isotopes are dependencies for do_react
    D_react = [vars_dep_res; vars_dep_react; vars_prop_react]
    PB.add_method_do!(
        rj,
        do_react,
        (
            PB.VarList_namedtuple_fields(fluxOceanBurial),
            PB.VarList_namedtuple(S_react), 
            PB.VarList_namedtuple(D_react),
        ),
        p=(CIsotopeType, SIsotopeType)
    )

    return rj
end


# Calculate rates
function do_react(m::PB.ReactionMethod, pars, (fluxOceanBurial, S, D), cellrange::PB.AbstractCellRange, deltat)

    (CIsotopeType, SIsotopeType) = m.p

    #%%%%% Reduced C species burial
    # Marine organic carbon burial
    if      pars.f_mocb[] == "original"
        D.mocb[]                = pars.k2_mocb[] * D.newp_relative[]^pars.f_mocb_b[]
    elseif  pars.f_mocb[] == "Uforced"
        D.mocb[]                = pars.k2_mocb[] * D.newp_relative[]^pars.f_mocb_b[] * D.UPLIFT[]
    elseif  pars.f_mocb[] == "O2dep"
        D.mocb[]                = pars.k2_mocb[] * D.newp_relative[]^pars.f_mocb_b[] * 2.1276*exp(-0.755*D.pO2PAL[])
    elseif  pars.f_mocb[] == "both"
        D.mocb[]                = pars.k2_mocb[] * D.newp_relative[]^pars.f_mocb_b[] * D.UPLIFT[] * 2.1276*exp(-0.755*D.pO2PAL[])
    else
        error("unknown f_mocb ", pars.f_mocb[])
    end
    if CIsotopeType <: PB.AbstractIsotopeScalar
        # delta of marine organic carbon burial
        D.mocb_delta[] = D.DIC_delta[] + D.D_mccb_DIC[] - D.D_B_mccb_mocb[]
    end
    mocb_isotope            = @PB.isotope_totaldelta(CIsotopeType, D.mocb[], D.mocb_delta[])
    
    # update tendencies and external fluxes
    S.O_sms[]               += D.mocb[]

    S.DIC_sms[]             += -mocb_isotope
    fluxOceanBurial.Corg[]  += mocb_isotope    

    #%%%% CP ratio
    if   pars.f_CPsea[] == "Fixed"
        D.CPsea[] = pars.CPsea0[]
    elseif pars.f_CPsea[] == "VCI"  # NB typo in Bergman (2004) has dependency reversed
        D.CPsea[] = (pars.f_CPsea_VCI_oxic[]*pars.f_CPsea_VCI_anoxic[] /
            ((1.0-D.ANOX[])*pars.f_CPsea_VCI_anoxic[] + D.ANOX[]*pars.f_CPsea_VCI_oxic[]))
    else
        error("unrecognized f_CPsea ", pars.f_CPsea[])
    end

    # Marine organic P burial
    mopb                    = (D.mocb[]/D.CPsea[])
    # Marine carbonate-associated P burial
    if      pars.f_capb[] == "original"
        capb                = pars.k7_capb[] * (D.newp_relative[]^pars.f_mocb_b[])
    elseif  pars.f_capb[] == "redox"
        capb                = (pars.k7_capb[] * (D.newp_relative[]^pars.f_mocb_b[])
                                    *(0.5+0.5*D.OX_relative[]))
    else
        error("unknown f_capb ", pars.f_capb[])
    end
    # Marine Fe-sorbed P burial NB: COPSE 5_14 uses copse_crash to limit at low P
    if      pars.f_fepb[] == "original"
        fepb = (pars.k6_fepb[]*D.OX_relative[] *
                PALEOcopse.COPSE.copse_crash(D.P_norm[], "fepb", D.tforce[]) )
    elseif  pars.f_fepb[] == "Dforced"
        fepb = D.DEGASS[]*(pars.k6_fepb[]*D.OX_relative[] *
                PALEOcopse.COPSE.copse_crash(D.P_norm[], "fepb", D.tforce[]) )
    elseif  pars.f_fepb[] == "sfw"
        fepb = D.sfw_relative[]*(pars.k6_fepb[]*D.OX_relative[] *
                                PALEOcopse.COPSE.copse_crash(D.P_norm[], "fepb", D.tforce[]) )
    elseif  pars.f_fepb[] == "pdep"
        fepb = pars.k6_fepb[]*D.OX_relative[]*D.P_norm[]
    else
        error("unknown f_fepb ", pars.f_fepb[])
    end

    totpb                   = mopb + capb + fepb    

    # update tendencies and external fluxes
    S.P_sms[] += -totpb
    fluxOceanBurial.Porg[]  += mopb
    fluxOceanBurial.Pauth[] += capb
    fluxOceanBurial.PFe[]   +=  fepb
    fluxOceanBurial.P[]     +=  totpb

    if pars.f_ncycle[]
        # Marine organic nitrogen burial
        D.monb[]                = D.mocb[]/pars.CNsea0[]
        # update tendencies, no external flux
        S.N_sms[]               += - D.monb[]
    end

    # Marine sulphur burial
    if pars.enableS[]
        if SIsotopeType <: PB.AbstractIsotopeScalar
            # Pyrite sulphur isotope fractionation relative to sulphate and gypsum
            if pars.f_sisotopefrac[] == "fixed"
                D.D_mpsb[] = 35.0
            elseif pars.f_sisotopefrac[] == "copse_O2"
                D.D_mpsb[] = 35.0*D.O_norm[]
            else
                error("unknown f_sisotopefrac ", pars.f_sisotopefrac[])
            end
        end

        # Marine gypsum sulphur burial
        mgsb                    = pars.k_mgsb[] * D.S_norm[] * D.CAL_norm[]
        mgsb_isotope            = @PB.isotope_totaldelta(SIsotopeType, mgsb, D.S_delta[])    

        # Marine pyrite sulphur burial
        if pars.f_pyrburial[] == "copse_noO2"  # dependent on sulphate and marine carbon burial
            mpsb                    = pars.k_mpsb[]*D.S_norm[]*(D.mocb[]/pars.k2_mocb[])
        elseif pars.f_pyrburial[] == "copse_O2"   # dependent on oxygen, sulphate, and marine carbon burial
            mpsb                    = pars.k_mpsb[]*D.S_norm[]/D.O_norm[]*(D.mocb[]/pars.k2_mocb[])
        else
            error("unknown f_pyrburial ", pars.f_pyrburial[])
        end
        mpsb_isotope            = @PB.isotope_totaldelta(SIsotopeType, mpsb, D.S_delta[] - D.D_mpsb[])

        # Update tendencies and external fluxes
        S.S_sms[]               += -(mgsb_isotope + mpsb_isotope)
        
        PB.add_if_available(S.CAL_sms, -mgsb)   # CAL is optional (eg not used in COPSE reloaded configs)
        S.O_sms[]               += 2*mpsb

        if pars.SRedoxAlk[]
            S.TAlk_sms[] += 2*D.mpsb[]
        end

        fluxOceanBurial.GYP[]   += mgsb_isotope
        fluxOceanBurial.PYR[]   += mpsb_isotope
    end
    
    return nothing
end


end # module
