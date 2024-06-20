# -*- coding: utf-8 -*-

"""
COPSE (Bergman, 2004) Phanerozoic biogeochemical model.
"""
module ModelBergman2004

import PALEOboxes as PB
using PALEOboxes.DocStrings

import PALEOcopse

"""
    ReactionModelBergman2004

Monolithic COPSE Bergman(2004) model, providing all biogeochemistry as a single Reaction in a single `global` Domain.

Requires reservoirs, forcings, and global temperature (ReactionGlobalTemperatureCK1992) to form a complete model.

# Parameters
$(PARS)

# Methods and Variables
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionModelBergman2004{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        ##### model constants (Table 4 Bergman etal 2004, consts.hh in COPSE 5_14 C code)
        PB.ParDouble("k1_oxfrac",  0.86,       units="",       description="initial oxic fraction"),
        PB.ParDouble("k2_mocb",    4.5e12,     units="mol/yr", description="ocean organic carbon burial"),
        PB.ParDouble("k3_nfix",    8.75e12,    units="mol/yr", description="nitrogen fixation"),
        PB.ParDouble("k4_denit",   4.3e12,     units="mol/yr", description="denitrification"),
        PB.ParDouble("k5_locb",        4.5e12,     units="mol/yr", description="land organic carbon burial"),
        PB.ParDouble("k6_fepb",    6e9,        units="mol/yr", description="Fe-P burial"),
        PB.ParDouble("k7_capb",    1.5e10,     units="mol/yr", description="Ca-P burial"),
        PB.ParDouble("k11_landfrac",   0.10345,    units="",       description="fract of weath P buried on land"),
        PB.ParDouble("k12_ccdeg",  6.65e12,    units="mol/yr",     description="carbonate C degassing"),
        PB.ParDouble("k13_ocdeg",  1.25e12,    units="mol/yr",     description="org C degassing"),
        PB.ParDouble("k14_carbw",      13.35e12,   units="mol/yr", description="carbonate weathering"),
        PB.ParDouble("k15_plantenhance",0.15,      units="",       description="weathering enhancement factor prior to vascular plant colonisation"),
        PB.ParDouble("k_fire",         100.0,      units="",       description="fire feedback"),
        PB.ParDouble("k16_PANtoO",     3.762,      units="",       description="For calculating mixing ratio O2 from normalised O2"),

  
        PB.ParDouble("k_silw",         NaN,        units="mol/yr", description="silicate weathering"),
        PB.ParDouble("k10_phosw",      NaN,        units="mol/yr", description="phosphorus weathering"),
        PB.ParDouble("k17_oxidw",      NaN,        units="mol/yr", description="oxidative org carbon weathering"),

        #### COPSE S system
        PB.ParDouble("k21_pyrw",       0.53e12,    units="mol S/yr",description="pyrite weathering"),
        PB.ParDouble("k22_gypw",       1e12,       units="mol S/yr",description="gypsum weathering"),       
        PB.ParDouble("k_pyrdeg",   0.0,        units="mol S/yr",   description="pyrite degassing"),
        PB.ParDouble("k_gypdeg",   0.0,        units="mol S/yr",   description="gypsum degassing"),
        PB.ParDouble("k_mpsb",     0.53e12,    units="mol S/yr",description="pyrite burial"),
        PB.ParDouble("k_mgsb",     1e12,       units="mol S/yr",description="gypsum burial"),


        ##########################################################################
        ####### Options controlling functional forms
        ##########################################################################

        # Organic carbon degassing
        # COPSE 5_14 C code (and Berman 2004) use "O2copsecrashprevent" which rolls
        # off organic carbon degassing at low pO2. This has a big effect at low pO2 when
        # oxidative weathering is oxygen-independent (ie Ordovician and earlier)
        PB.ParString("f_ocdeg",    "O2indep",  allowed_values=["O2indep", "O2copsecrashprevent"], description="roll of orgC degassing at low pO2"),
        PB.ParString("f_temp", "CK1992",   description="functional form of temperature estimate"),

        # C (carbonate reservoir) dependence of carbonate weathering
        PB.ParString("f_carbwC",       "Cindep",   allowed_values=["Cindep", "Cprop"], description="C (carbonate reservoir) dependence of carbonate weathering"),

        ####### Oxidative weathering
        PB.ParString("f_oxwO",         "PowerO2",  allowed_values=["PowerO2", "SatO2"], description="oxidative weathering functional form"),
        PB.ParDouble("f_oxw_a",        0.5,        units="", description="oxidative weathering dependency on O2 concentration"),
        #f_oxw_halfsat  : 1e-6 #set small value to represent O2-indep with limit to satisfy ODE integrator

        # Marine N cycle
        PB.ParDouble("f_nfix_power", 2.0,      units="",       description="SD update to BM Matlab code which has f_nfix_power : 1. This fixes discrepancy between P (phosphorus) and COPSE 5_14 (most visible in Cambrian)"),
        PB.ParString("f_nfix_nreplete","Off",  allowed_values=["Off", "Sign"], description="Bug compatitibility with COPSE 5_14 C code, which has Sign)"),

        # Marine ecosystem
        PB.ParDouble("newp0",      225.956),

        # Marine CPN ratio
        PB.ParString("f_CPsea",    "Fixed",    allowed_values=["Fixed", "VCI"], description="Functional form of CPsea ratio"),
        PB.ParDouble("CPsea0",      250.0,     description="for f_CPsea:\"Fixed\""),
        # Van Cappellen & Ingall (f_CPsea:"VCI") marine C/P burial ratio consts
        # These are Redfield Revisited 1 steady-state values
        PB.ParDouble("f_CPsea_VCI_oxic",217.0,     description="Van Cappellen & Ingall (f_CPsea:\"VCI\") marine C/P burial ratio consts"),
        PB.ParDouble("f_CPsea_VCI_anoxic",4340.0,  description="Van Cappellen & Ingall (f_CPsea:\"VCI\") marine C/P burial ratio consts"),
        #f_CPsea_VCI_oxic : 200.0     # Original VC&I values
        #f_CPsea_VCI_anoxic : 4000.0

        PB.ParDouble("CNsea0",      37.5,       description="Always fixed"),

        PB.ParDouble("f_mocb_b",    2.0,        description="marine organic carbon burial power-law dependency on new production"),

        #Marine pyrite sulphur burial dependency on oxygen
        PB.ParString("f_pyrburial", "copse_O2", allowed_values=["copse_O2", "copse_noO2"], description="Marine pyrite sulphur burial dependency on oxygen"),

        # C isotope fractionation calculation
        PB.ParString("f_cisotopefrac", "copse_base", allowed_values=["fixed", "copse_base", "copse_noO2"],
                                                description="fractionation calculation method"),
        # S isotope fractionation calculation
        PB.ParString("f_sisotopefrac","fixed",  allowed_values=["fixed", "copse_O2"], description="S isotope fractionation calculation.")

    )

end


function PB.register_methods!(rj::ReactionModelBergman2004)
 
    vars_res, vars_sms, vars_dep_res     = PB.Reservoirs.ReservoirLinksVector(
        Dict("CIsotope"=>PB.IsotopeLinear, "SIsotope"=>PB.IsotopeLinear),
        [
            ("P",           "mol",    "Marine phosphorus"),
            ("N",           "mol",    "Marine nitrogen"),
            ("O",           "mol",   "Atm-ocean oxygen"),
            ("C::CIsotope",  "mol",    "Sedimentary carbonate"),
            ("G::CIsotope",  "mol",    "Sedimentary organic carbon"),
            ("A::CIsotope",  "mol",    "Atm-ocean inorganic carbon"),
            ("CAL",         "mol",   "Marine calcium"),
            ("GYP::SIsotope","mol",    "Sedimentary gypsum"),
            ("PYR::SIsotope","mol",    "Sedimentary pyrite"),
            ("S::SIsotope",  "mol",    "Marine sulphate")
        ])
 
    

    # COPSE forcings - dependencies, assumed required in do_stateandeqb
    vars_dep_stateandeqb = PB.VarVector(PB.VarDepScalar, 
        [
            ("tforce",  "yr"   , "time for external forcings"),
            ("SOLAR",   "W m-2", "solar luminosity"),
            ("UPLIFT",  "",      "uplift scaling"),
            ("DEGASS",  "",      "degassing scaling"),
            ("W",       "",      "plant weathering scaling"),
            ("EVO",     "",      "plant evolution scaling"),
            ("Bforcing","",      "calcerous plankton evolution"),
            ("CPland_relative","","land CP burial ratio scaling"),
            # include TEMP as that is supplied by external Reaction
            ("TEMP",    "K",      "global temperature")
        ])
    # Properties we calculate in do_stateandeqb
    vars_prop_stateandeqb = PB.VarVector(PB.VarPropScalar,
        [
            ("pO2PAL",  "",           "atmospheric oxygen relative to present-day"),
            ("pCO2PAL", "",           "atmospheric pCO2 relative to present-day"),
            ("pCO2atm", "atm",        "atmospheric pCO2"),
            ("d_mccb", "per mil",     "D13C fractionation carbonate - atm-ocean inorganic carbon"),
            ("d_locb", "per mil",     "D13C fractionation land organic carbon burial - atm-ocean inorganic carbon"),
            ("d_mocb", "per mil",     "D13C fractionation marine organic carbon burial - atm-ocean inorganic carbon"),
            ("mccb_delta", "per mil", "D13C fractionation of marine calcium carbonate burial"),
              # sulphur isotopes
            ("D_mpsb", "per mil",     "D34S fractionation pyrite - marine sulphate")
        ])

    # Properties we calculate in do_react
    vars_prop_react = PB.VarVector(PB.VarPropScalar,
        [
            # copse_landbiota
            ("V_T",     "",           "effect of temp on VEG"),
            ("V_co2",   "",           "effect of CO2 on VEG"),
            ("V_o2",    "",           "effect of O2 on VEG"),
            ("V_npp",   "",           "full VEG limitation"),
            ("mrO2",    "",           "atmospheric O2 mixing ratio"),
            ("ignit",   "",           "ignition probability"),
            ("firef",   "",           "effect of fire on VEG"),
            ("VEG",     "",           "Mass of terrestrial biosphere"),
            # weathering
            ("f_preplant","",         "Pre-plant silicate weathering T and pCO2 factor"),
            ("f_plant", "",           "Plant silicate weathering T and pCO2 factor"),
            ("g_preplant","",         "Pre-plant carbonate weathering T and pCO2 factor"),
            ("g_plant", "",           "Plant carbonate weathering T and pCO2 factor"),
            ("VWmin",   "",           "Plant weighting for weathering T factors"),
            ("f_co2",   "",           "Plant weighted silicate weathering T and pCO2 factor"),
            ("g_co2",   "",           "Plant weighted carbonate weathering T and pCO2 factor"),
            ("w_plantenhance","",      "Plant weighted weathering rate factor"),
            ("silw",    "molC/yr",    "Silicate weathering"),
            ("carbw",   "molC/yr",    "Carbonate weathering"),
            ("oxw_fac", "",           "Oxidative weathering oxygen-dependent factor"),
            ("oxidw",   "molC/yr",    "Oxidative organic C weathering"),
            ("gypw",    "molS/yr",    "Gypsum weathering"),
            ("pyrw",    "molS/yr",    "Pyrite weathering"),
            ("phosw_s", "molP/yr",    "Phosphorus weathering, silicate-associated"),
            ("phosw_c", "molP/yr",    "Phosphorus weathering, carbonate-associated"),
            ("phosw_o", "molP/yr",    "Phosphorus weathering, organic carbon-associated"),
            ("phosw",   "molP/yr",    "Phosphorus weathering, total"),
            ("pland",   "molP/yr",    "Phosphorus weathering to land"),
            ("psea",    "molP/yr",    "Phosphorus weathering to sea"),
            # degassing
            ("ccdeg",   "molC/yr",    "Carbonate degassing"),
            ("ocdeg",   "molC/yr",    "Organic carbon degassing"),
            ("gypdeg",  "molS/yr",    "Gypsum degassing"),
            ("pyrdeg",  "molS/yr",    "Pyrite degassing"),
            # copse_marinebiota
            ("Pconc",   "molP/kgsw",  "marine P conc"),
            ("Nconc",   "molN/kgsw",  "marine N conc"),
            ("newp",    "",           "marine new production"),
            ("ANOX",    "",           "marine anoxic fraction"),
            ("CPsea",   "",           "marine C:P ratio"),
            ("nfix",    "molN/yr",    "marine nitrogen fixation"),
            ("denit",   "molP/yr",    "marine denitrification"),
            # burial
            ("mccb",   "molC/yr",     "Carbonate burial"),
            ("mocb",   "molC/yr",     "Marine organic carbon burial"),
            ("locb",   "molC/yr",     "Land organic carbon burial"),
            ("mopb",   "molP/yr",     "Marine organic P burial"),
            ("monb",   "molN/yr",     "Marine organic N burial"),
            ("capb",   "molP/yr",     "Marine carbonate-associated P burial"),
            ("fepb",   "molP/yr",     "Marine Fe-sorbed P burial"),
            ("mgsb",   "molS/yr",     "Marine gypsum sulphur burial"),
            ("mpsb",   "molS/yr",     "Marine pyrite sulphur burial")
        ])
    
    # Now assemble these lists into properties and dependencies for each method,
    # following the Matlab COPSE convention that 
    #   S - reservoir (a VarDep) and reservoir_sms (a VarContrib)
    #   D - dependencies (VarDep) and properties to calculate (VarProp)

    S_stateandeqb = vars_res
    D_stateandeqb = [vars_dep_res; vars_dep_stateandeqb; vars_prop_stateandeqb]
    PB.add_method_do!(
        rj,
        do_stateandeqb,
        (PB.VarList_namedtuple(S_stateandeqb), PB.VarList_namedtuple(D_stateandeqb)),
    )

    S_react = vars_res
    D_react = [vars_dep_res; vars_dep_stateandeqb; PB.VarDep.(vars_prop_stateandeqb); vars_prop_react]
    PB.add_method_do!(
        rj,
        do_react,
        (PB.VarList_namedtuple(S_react), PB.VarList_namedtuple(D_react)),
    )

    S_updatestate = [vars_res; vars_sms]
    D_updatestate = [vars_dep_res; vars_dep_stateandeqb; PB.VarDep.(vars_prop_stateandeqb); PB.VarDep.(vars_prop_react)]
    PB.add_method_do!(
        rj,
        do_updatestate,
        (PB.VarList_namedtuple(S_updatestate), PB.VarList_namedtuple(D_updatestate)),
    )

    return nothing
end

function set_parameters_modern_steady_state(copsemodel::ReactionModelBergman2004)
    # set parameters for steady-state
    pars = copsemodel.pars
    PB.setvalue!(pars.k17_oxidw, pars.k2_mocb[] + pars.k5_locb[] - pars.k13_ocdeg[])
    @info "set" pars.k17_oxidw

    PB.setvalue!(pars.k_silw, -pars.k2_mocb[] - pars.k5_locb[]
                    + pars.k17_oxidw[] + pars.k13_ocdeg[] + pars.k12_ccdeg[])
    @info "set" pars.k_silw

    PB.setvalue!(pars.k10_phosw, (pars.k2_mocb[]/pars.CPsea0[]  + pars.k7_capb[] + pars.k6_fepb[]) /
                            (1.0-pars.k11_landfrac[]))
    @info "set" pars.k10_phosw

end




function do_stateandeqb(m::PB.ReactionMethod, pars, (S, D), cellrange::PB.AbstractCellRange, deltat)

    # Add user-friendly atmospheric concentrations etc
    D.pO2PAL[] = D.O_norm[]
    D.pCO2PAL[] = D.A_norm[]  # pre-industrial = 1
    D.pCO2atm[] = D.pCO2PAL[]*PB.Constants.k_preindCO2atm   # pre-industrial = 280e-6


    # Carbon isotope fractionation (relative to total CO2 (A reservoir)
    if pars.f_cisotopefrac[] == "fixed"
        D.d_mocb[] = -30.0
        D.d_locb[] = -30.0
        D.d_mccb[] = 0.0
    elseif pars.f_cisotopefrac[] == "copse_base"
        D.d_locb[], D_P, D.d_mocb[], D_B, D.d_mccb[], d_ocean, d_atmos =
             copse_Cisotopefrac( D.TEMP[], D.pCO2PAL[], D.O_norm[] )
    elseif pars.f_cisotopefrac[] == "copse_noO2"
        D.d_locb[], D_P, D.d_mocb[], D_B, D.d_mccb[], d_ocean, d_atmos =
             copse_Cisotopefrac( D.TEMP[], D.pCO2PAL[], 1.0 )
    else
        error("unknown f_cisotopefrac ", pars.f_cisotopefrac[])
    end
    # delta of marine carbonate burial
    D.mccb_delta[] = D.A_delta[] + D.d_mccb[]

    # Pyrite sulphur isotope fractionation relative to sulphate and gypsum
    if   pars.f_sisotopefrac[] == "fixed"
        D.D_mpsb[] = 35.0
    elseif pars.f_sisotopefrac[] == "copse_O2"
        D.D_mpsb[] = 35.0*D.O_norm[]
    else
        error("unknown f_sisotopefrac ", pars.f_sisotopefrac[])
    end

    return nothing
end


function do_react(m::PB.ReactionMethod, pars, (S, D), cellrange::PB.AbstractCellRange, deltat)

    # land biota
    copse_landbiota_bergman2004( pars, S, D, D.TEMP[] )

    #%%%%%%%% calculate weathering
    #%%%% Plant effects on weathering

    D.f_preplant[] = copse_f_T(D.TEMP[]) * max(D.pCO2PAL[], 1e-3)^0.5
    D.f_plant[]    = copse_f_T(D.TEMP[]) * ( 2.0*max(D.pCO2PAL[], 1e-3) / (1.0 + max(D.pCO2PAL[], 1e-3)) )^0.4

    D.g_preplant[] = copse_g_T(D.TEMP[]) * max(D.pCO2PAL[], 1e-3)^0.5
    D.g_plant[]    = copse_g_T(D.TEMP[]) * ( 2.0*max(D.pCO2PAL[], 1e-3) / (1.0 + max(D.pCO2PAL[], 1e-3)) )^0.4

    D.VWmin[] = min(D.VEG[]*D.W[], 1.0)

    D.f_co2[] = D.f_preplant[]*(1.0 - D.VWmin[]) + D.f_plant[]*D.VWmin[]
    D.g_co2[] = D.g_preplant[]*(1.0 - D.VWmin[]) + D.g_plant[]*D.VWmin[]

    D.w_plantenhance[] = (pars.k15_plantenhance[] +
                (1.0 - pars.k15_plantenhance[]) * D.W[] * D.VEG[])

    # silicate and carbonate weathering
    D.silw[] = pars.k_silw[]*D.UPLIFT[]*D.w_plantenhance[] * D.f_co2[]

    carbw_fac = D.UPLIFT[]*D.w_plantenhance[]*D.g_co2[]
    if pars.f_carbwC[] == "Cindep"   # Copse 5_14
        D.carbw[] = pars.k14_carbw[]*carbw_fac
    elseif pars.f_carbwC[] == "Cprop"   # A generalization for varying-size C reservoir
        D.carbw[] = pars.k14_carbw[]*carbw_fac*D.C_norm[]
    else
       error("Unknown f_carbwC ", pars.f_carbwC[])
    end

    #%%%% Oxidative weathering
    # Functional form of oxidative weathering
    if pars.f_oxwO[] == "PowerO2"   # Copse 5_14 base with f_oxw_a = 0.5
        D.oxw_fac[] = D.O_norm[]^pars.f_oxw_a[]
    elseif pars.f_oxwO[] == "SatO2"
        D.oxw_fac[] = D.O_norm[] /(D.O_norm[] + pars.f_oxw_halfsat[])
    else
       error("Unknown f_foxwO ", pars.f_oxwO[])
    end
    # C oxidative weathering
    D.oxidw[] = pars.k17_oxidw[] * D.UPLIFT[]*D.G_norm[]*D.oxw_fac[]

    # Sulphur weathering

    # Gypsum weathering tied to carbonate weathering
    D.gypw[] = pars.k22_gypw[] * D.GYP_norm[] * carbw_fac
    # Pyrite oxidative weathering with same functional form as carbon
    D.pyrw[] = pars.k21_pyrw[] * D.UPLIFT[]*D.PYR_norm[] * D.oxw_fac[]

    # P weathering and delivery to land and sea
    D.phosw_s[] = pars.k10_phosw[] * (2.0/12.0)*(D.silw[]/pars.k_silw[])
    D.phosw_c[] = pars.k10_phosw[] * (5.0/12.0)*(D.carbw[]/pars.k14_carbw[])
    D.phosw_o[] = pars.k10_phosw[] * (5.0/12.0)*(D.oxidw[]/pars.k17_oxidw[])
    D.phosw[]   = D.phosw_s[] + D.phosw_c[] + D.phosw_o[]

    D.pland[]   = pars.k11_landfrac[]*D.VEG[]*D.phosw[]
    pland0      = pars.k11_landfrac[]*pars.k10_phosw[]

    D.psea[]    = D.phosw[] - D.pland[]

    #%%%%%%%% calculate degassing
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    #  Inorganic carbon
    D.ccdeg[]  = pars.k12_ccdeg[]*D.DEGASS[]*D.C_norm[]*D.Bforcing[]

    # Organic carbon
    ocdeg_raw = pars.k13_ocdeg[]*D.DEGASS[]*D.G_norm[]
    if pars.f_ocdeg[] == "O2indep"
        D.ocdeg[] = ocdeg_raw
    elseif pars.f_ocdeg[] == "O2copsecrashprevent"
        # COPSE 5_14 does this (always) apparently to prevent pO2 dropping to zero ?
        # This has a big effect when pO2 dependence of oxidative weathering switched off
        D.ocdeg[] = ocdeg_raw*PALEOcopse.COPSE.copse_crash(D.O_norm[], "ocdeg", D.tforce[])
    else
        error("unrecogized pars.f_ocdeg ", pars.f_ocdeg[])
    end

    # Sulphur
    D.pyrdeg[] = pars.k_pyrdeg[]*D.PYR_norm[]*D.DEGASS[]
    D.gypdeg[] = pars.k_gypdeg[]*D.GYP_norm[]*D.DEGASS[]


    #%%%%%%% Marine biota
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    copse_marinebiota_bergman2004(pars, D.tforce[], S, D )

    #%%%%%%% Burial
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    #%%%%% Oxidised C species burial
    D.mccb[] = D.carbw[] + D.silw[]  # disguised alkalinity balance

    #%%%%% Reduced C species burial
    # Marine organic carbon burial
    D.mocb[] = pars.k2_mocb[] * (D.newp[]/pars.newp0[])^pars.f_mocb_b[]
    # Land organic carbon burial
    D.locb[] = pars.k5_locb[] * (D.pland[]/pland0) * D.CPland_relative[]

    # Marine organic P burial
    D.mopb[] = (D.mocb[]/D.CPsea[])
    # Marine carbonate-associated P burial
    D.capb[] = pars.k7_capb[] * ((D.newp[]/pars.newp0[])^pars.f_mocb_b[])
    # Marine Fe-sorbed P burial NB: COPSE 5_14 uses PALEOcopse.COPSE.copse_crash to limit at low P
    D.fepb[] = ((pars.k6_fepb[]/pars.k1_oxfrac[]) * (1.0 - D.ANOX[]) *
                    PALEOcopse.COPSE.copse_crash(D.P_norm[], "fepb", D.tforce[]) )

    # Marine organic nitrogen burial
    D.monb[] = D.mocb[]/pars.CNsea0[]

    # Marine sulphur burial
    # Marine gypsum sulphur burial
    D.mgsb[] = pars.k_mgsb[] * D.S_norm[] * D.CAL_norm[]
    # Marine pyrite sulphur burial
    if pars.f_pyrburial[] == "copse_noO2"  # dependent on sulphate and marine carbon burial
        D.mpsb[] = pars.k_mpsb[]*D.S_norm[]*(D.mocb[]/pars.k2_mocb[])
    elseif pars.f_pyrburial[] == "copse_O2"   # dependent on oxygen, sulphate, and marine carbon burial
        D.mpsb[] = pars.k_mpsb[]*D.S_norm[]/D.O_norm[]*(D.mocb[]/pars.k2_mocb[])
    else
        error("unknown f_pyrburial ", pars.f_pyrburial[])
    end

    return nothing
end


function do_updatestate(m::PB.ReactionMethod, (S, D), cellrange::PB.AbstractCellRange, deltat)

    # %%%% Atmosphere / ocean reservoirs
     # Oxygen
    S.O_sms[] = (D.locb[] + D.mocb[] - D.oxidw[] - D.ocdeg[] +
                        2.0*(D.mpsb[] - D.pyrw[] - D.pyrdeg[]) )

    # Carbon
    S.A_sms.v[] = (-D.locb[] - D.mocb[] + D.oxidw[] + D.ocdeg[] + D.ccdeg[] +
                    D.carbw[] - D.mccb[])

    # Marine nutrient reserviors
    S.P_sms[] = D.psea[] - D.mopb[] - D.capb[] - D.fepb[]
    S.N_sms[] = D.nfix[] - D.denit[] - D.monb[]

    # Marine calcium
    S.CAL_sms[] = D.silw[] + D.carbw[] + D.gypw[] - D.mccb[] - D.mgsb[]

    # Crustal C reservoirs
    # Buried organic C
    S.G_sms.v[] = D.locb[] + D.mocb[] - D.oxidw[] - D.ocdeg[]

    # Buried carb C
    S.C_sms.v[] = D.mccb[] - D.carbw[] - D.ccdeg[]

    #%%%%%%% calculate isotopic fractionation of reservoirs
    # deltaORG_C*ORG_C
    S.G_sms.v_moldelta[] =  (D.mocb[]*( D.A_delta[] + D.d_mocb[] )
                        + D.locb[]*( D.A_delta[] + D.d_locb[] )
                        - D.oxidw[]*D.G_delta[]  - D.ocdeg[]*D.G_delta[])

    # deltaCARB_C*CARB_C
    S.C_sms.v_moldelta[] =  (D.mccb[]*D.mccb_delta[]
                           -  D.carbw[]*D.C_delta[] - D.ccdeg[]*D.C_delta[] )

    # delta_A * A
    S.A_sms.v_moldelta[] = (-D.locb[]*( D.A_delta[] + D.d_locb[] )
                        - D.mocb[]*( D.A_delta[] + D.d_mocb[] )
                        + D.oxidw[]*D.G_delta[] + D.ocdeg[]*D.G_delta[]
                        + D.ccdeg[]*D.C_delta[] + D.carbw[]*D.C_delta[]
                        - D.mccb[]*(D.A_delta[] + D.d_mccb[]))

    # Sulphur

    # Marine sulphate
    S.S_sms.v[] = D.gypw[] + D.pyrw[] - D.mgsb[] - D.mpsb[] +D.gypdeg[] + D.pyrdeg[]
    # Buried pyrite S
    S.PYR_sms.v[] = D.mpsb[] - D.pyrw[] - D.pyrdeg[]
    # Buried gypsum S
    S.GYP_sms.v[] = D.mgsb[] - D.gypw[] -D.gypdeg[]

    # Isotopes

    # deltaPYR_S*PYR_S
    S.PYR_sms.v_moldelta[] =  (D.mpsb[]*(D.S_delta[] - D.D_mpsb[])
                            - D.pyrw[]*D.PYR_delta[] - D.pyrdeg[]*D.PYR_delta[])
    # deltaGYP_S*GYP_S
    S.GYP_sms.v_moldelta[] =  (D.mgsb[]*D.S_delta[]
                            - D.gypw[]*D.GYP_delta[] - D.gypdeg[]*D.GYP_delta[])
    # delta_S * S
    S.S_sms.v_moldelta[] =    (D.gypw[]*D.GYP_delta[] + D.pyrw[]*D.PYR_delta[]
                            - D.mgsb[]*D.S_delta[] - D.mpsb[]*(D.S_delta[] - D.D_mpsb[] )
                            + D.gypdeg[]*D.GYP_delta[] + D.pyrdeg[]*D.PYR_delta[])

    return nothing
end


"""
    copse_landbiota_bergman2004(pars, S, D, TEMP)

COPSE OCT dynamic land vegetation model
"""
function copse_landbiota_bergman2004(pars, S, D, TEMP)
    

    # effect of temp on VEG
    # SD update to BM Matlab - this fixes absence of T limitation
    # visible as slightly high VEG
    D.V_T[] = 1.0 - (((TEMP - PB.Constants.k_CtoK-25.0)/25.0)^2)

    # effect of CO2 on VEG
    P_atm = D.pCO2atm[]*1.0e6
    P_half = 183.6
    P_min = 10.0
    D.V_co2[] = (P_atm - P_min) / (P_half + P_atm - P_min)

    # effect of O2 on VEG
    # SD update to BM Matlab - this fixes an incorrect use of mrO2
    # visible as slightly high VEG
    D.V_o2[] = 1.5 - 0.5*D.pO2PAL[]

    # full VEG limitation
    D.V_npp[] = 2.0*D.EVO[]*D.V_T[]*D.V_o2[]*D.V_co2[]

    # fire feedback
    # calculate atmospheric mixing ratio of O2 (for constant atmospheric N etc!)
    # (only used for fire ignition probability)
    D.mrO2[] = D.pO2PAL[] / (D.pO2PAL[] + PB.Constants.k16_PANtoO)
    D.ignit[] = max(586.2*D.mrO2[] - 122.102, 0.0)
    D.firef[] = pars.k_fire[]/(pars.k_fire[] - 1.0 + D.ignit[])

    # Mass of terrestrial biosphere
    D.VEG[] = D.V_npp[] * D.firef[]
    return nothing
end

"""
    copse_marinebiota_bergman2004(pars, tmodel, S, D )

COPSE marine ecosystem model    
"""
function copse_marinebiota_bergman2004(pars, tmodel, S, D )
    

    # convert marine nutrient reservoir moles to micromoles/kg concentration
    D.Pconc[] = D.P_norm[] * 2.2
    D.Nconc[] = D.N_norm[] * 30.9 ;
    #  clunky way to get normalization values
    # fails as Reaction doesn't see norm_value for variable it is linked to
    # P0 = S.P.norm_value
    # N0 = S.N.norm_value
    P0 = S.P[]/D.P_norm[]
    N0 = S.N[]/D.N_norm[]

    #%%%%% Marine new production
    D.newp[] = 117.0 * min(D.Nconc[]/16.0, D.Pconc[])

    #%%%%%% OCEAN ANOXIC FRACTION

    D.ANOX[] = max( 1.0 - pars.k1_oxfrac[]*D.pO2PAL[] * pars.newp0[]/D.newp[], 0.0 )

    #%%%% CP ratio
    if   pars.f_CPsea[] == "Fixed"
        D.CPsea[] = pars.CPsea0[]
    elseif pars.f_CPsea[] == "VCI"  # NB typo in Bergman (2004) has dependency reversed
        D.CPsea[] = (pars.f_CPsea_VCI_oxic[]*pars.f_CPsea_VCI_anoxic[] /
         ((1.0-D.ANOX[])*pars.f_CPsea_VCI_anoxic[] + D.ANOX[]*pars.f_CPsea_VCI_oxic[]))
    else
        error("unrecognized f_CPsea ", pars.f_CPsea[])
    end

    #%%%%% nitrogen cycle
    if (S.N[]/16.0) < S.P[]
        D.nfix[] = pars.k3_nfix[] *( ( S.P[] - S.N[]/16.0 ) / ( P0  - N0/16.0) )^pars.f_nfix_power[]
    else
        if   pars.f_nfix_nreplete[] == "Off" # Surely more defensible ?
            D.nfix[] = 0.0
        elseif pars.f_nfix_nreplete[] == "Sign" # SD - COPSE 5_14 C code has this (?!)
            D.nfix[] = pars.k3_nfix[] *( - ( S.P[] - S.N[]/16.0 ) / ( P0  - N0/16.0)  )^pars.f_nfix_power[]
            print("COPSE_equations -ve nfix (check pars.f_nfix_nreplete) tmodel ", tmodel)
        else
            error("unrecognized f_nfix_nreplete ", pars.f_nfix_nreplete[])
        end
    end
    # Denitrification NB: COPSE 5_14 uses PALEOcopse.COPSE.copse_crash to limit at low N
    D.denit[] = (pars.k4_denit[] * (1.0 + D.ANOX[] / (1.0 - pars.k1_oxfrac[]) ) *
                     PALEOcopse.COPSE.copse_crash(D.N_norm[], "denit", tmodel))
    return nothing
end

""" 
    copse_Cisotopefrac(Tkelvin, pCO2PAL, pO2PAL [, phi])

Carbon isotope fractionation from Bergman(2004) COPSE model
"""
function copse_Cisotopefrac(Tkelvin, pCO2PAL, pO2PAL, phi=nothing)
    
    if isnothing(phi)
        phi = 0.01614  # fraction of C in atmosphere:ocean
    end

    # ocean total dissolved CO2
    d_ocean = phi*(9483.0/Tkelvin-23.89)
    # atmosphere CO2
    d_atmos = (phi-1.0)*(9483.0/Tkelvin-23.89)

    # marine calcite burial
    d_mccb = d_ocean + 15.10 - 4232.0/Tkelvin

    # fractionation between marine organic and calcite burial
    D_B = 33.0 - 9.0/sqrt(max(pCO2PAL, 1e-3)) + 5.0*(pO2PAL - 1.0)
    # marine organic carbon burial
    d_mocb = d_mccb - D_B

    # fractionation between terrestrial organic burial and atmospheric CO2
    D_P = 19.0 + 5.0*(pO2PAL-1.0)
    d_locb = d_atmos - D_P

    return (d_locb, D_P, d_mocb, D_B, d_mccb, d_ocean, d_atmos)
end

""" 
    copse_f_T(Tk)

COPSE_F_T Temperature dependence of silicate weathering from COPSE (Bergman 2004)
Compute `f_T` weathering rate (relative to Tzero = 15 C) at temperature `Tk` Kelvin
"""
function copse_f_T(Tk)
    

    Tzero = PB.Constants.k_CtoK+15.0

    f_T_arg = 1.0 + 0.038*(Tk - Tzero)
    if (f_T_arg < 0.0)
        print("copse_f_T: negative f_T_arg ", f_T_arg," for Tk ", Tk, " resetting to 0.0\n")
        f_T_arg = zero(Tk)
    end

    copse_f_T = exp(0.09*(Tk-Tzero)) * f_T_arg^0.65

    return copse_f_T
end

"""
    copse_g_T(Tk) 

Temperature dependence of carbonate weathering from COPSE (Bergman 2004)
Compute `g_T` carbonate weathering rate relative to Tzero = 15 C
NB: will go negative to T <~ 5C 

# Arguments
- `Tk`: Temperature (Kelvin)               
"""
function copse_g_T(Tk)
    

    Tzero = PB.Constants.k_CtoK+15.0

    g_T = 1.0 + 0.087*(Tk - Tzero)
    if (g_T < 0.0)
        print("copse_g_T: negative g_T ", g_T," for Tk ", Tk, " resetting to 0.0\n")
        g_T = zero(Tk)
    end

    return g_T
end



end # module
