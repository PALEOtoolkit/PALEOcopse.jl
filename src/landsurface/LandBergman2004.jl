# -*- coding: utf-8 -*-
module LandBergman2004


import PALEOboxes as PB
using PALEOboxes.DocStrings

"""
    ReactionLandBergman2004

COPSE Bergman(2004) land surface

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionLandBergman2004{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        ##### model constants (Table 4 Bergman etal 2004, consts.hh in COPSE 5_14 C code)
        PB.ParDouble("k5_locb",        4.5e12,     units="mol/yr", description="land organic carbon burial"),
        PB.ParDouble("k11_landfrac",   0.10345,    units="",       description="fract of weath P buried on land"),
        PB.ParDouble("k14_carbw",      13.35e12,   units="mol/yr", description="carbonate weathering"),
        PB.ParDouble("k15_plantenhance",0.15,      units="",       description="weathering enhancement factor prior to vascular plant colonisation"),
        PB.ParDouble("k_fire",         100.0,      units="",       description="fire feedback"),
        PB.ParDouble("k16_PANtoO",     3.762,      units="",       description="For calculating mixing ratio O2 from normalised O2"),

        PB.ParDouble("k_silw",         NaN,        units="mol/yr", description="silicate weathering"),
        PB.ParDouble("k10_phosw",      NaN,        units="mol/yr", description="phosphorus weathering"),
        PB.ParDouble("k17_oxidw",      NaN,        units="mol/yr", description="oxidative org carbon weathering"),

        #### COPSE S system
        PB.ParBool("enableS",          true,                        description="enable S weathering"),
        PB.ParDouble("k21_pyrw",       0.53e12,    units="mol S/yr",description="pyrite weathering"),
        PB.ParDouble("k22_gypw",       1e12,       units="mol S/yr",description="gypsum weathering"),       


        ##########################################################################
        ####### Options controlling functional forms
        ##########################################################################

        # C (carbonate reservoir) dependence of carbonate weathering
        PB.ParString("f_carbwC",       "Cindep",   allowed_values=["Cindep", "Cprop"], description="C (carbonate reservoir) dependence of carbonate weathering"),

        ####### Oxidative weathering
        PB.ParString("f_oxwO",         "PowerO2",  allowed_values=["PowerO2", "SatO2"], description="oxidative weathering functional form"),
        PB.ParDouble("f_oxw_a",        0.5,        units="", description="oxidative weathering dependency on O2 concentration"),
        #f_oxw_halfsat  : 1e-6 #set small value to represent O2-indep with limit to satisfy ODE integrator
        
        # land biota
        PB.ParString("f_landbiota", "Dynamic",  allowed_values=["Dynamic", "Prescribed"],
            description="enable/disable dynamic land biota (VEG required as forcing if Prescribed)"),

        # Organic carbon burial
        PB.ParString("f_locb",      "original", allowed_values=["original", "Prescribed"]),

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

function PB.register_methods!(rj::ReactionLandBergman2004)

    if rj.pars.f_landbiota[] == "Dynamic"
        var_VEG = PB.VarPropScalar("VEG",     "",           "Mass of terrestrial biosphere")
    elseif rj.pars.f_landbiota[] == "Prescribed"
        # VEG supplied as forcing
        var_VEG = PB.VarDepScalar("global.VEG",     "",           "Mass of terrestrial biosphere")
    else
        error("unknown f_landbiota='$(pars.f_landbiota[])'")
    end

    # dependencies
    vars_dep = PB.VarVector(PB.VarDepScalar,
        [
            # Forcings
            ("global.tforce",  "yr"   , "time for external forcings"),
            ("global.UPLIFT",  "",      "uplift scaling"),
            ("global.W",       "",      "plant weathering scaling"),
            ("(global.EVO)",     "",      "plant evolution scaling"),           
            ("(global.CPland_relative)","","land CP burial ratio scaling"),
            ("(global.locbpert)",    "",  "land organic carbon burial forcing (f_locb==Prescribed only)"),
            # 
            ("global.TEMP",    "K",      "global temperature"),
            #
            ("(D_P_CO2_locb)", "per mil",  "d13C fractionation between terrestrial organic burial and atmospheric CO2"),
            # atm
            ("atm.pO2PAL",  "",           "atmospheric oxygen relative to present-day"),
            ("atm.pCO2PAL", "",           "atmospheric pCO2 relative to present-day"),
            ("atm.pCO2atm", "ppm",        "atmospheric pCO2"),
            ("(atm.CO2_delta)", "per mil",  "atmospheric pCO2 delta 13C"),
            # sedcrust
            ("sedcrust.C_norm",        "",    "Sedimentary carbonate normalized to present day"),
            ("(sedcrust.C_delta)",        "",    "Sedimentary carbonate d13C"),
            ("sedcrust.G_norm",        "",    "Sedimentary organic carbon normalized to present day"),
            ("(sedcrust.G_delta)",        "",    "Sedimentary organic carbon d13C"),
            ("(sedcrust.GYP_norm)",      "",    "Sedimentary gypsum normalized to present day"),
            ("(sedcrust.GYP_delta)",      "",    "Sedimentary gypsum d34S"),
            ("(sedcrust.PYR_norm)",      "",    "Sedimentary pyrite normalized to present day"),
            ("(sedcrust.PYR_delta)",      "",    "Sedimentary pyrite d34S")
        ])

    # properties we calculate
    vars_prop = PB.VarVector(PB.VarPropScalar,
        [
            # copse_landbiota
            ("V_T",     "",           "effect of temp on VEG"),
            ("V_co2",   "",           "effect of CO2 on VEG"),
            ("V_o2",    "",           "effect of O2 on VEG"),
            ("V_npp",   "",           "full VEG limitation"),
            ("mrO2",    "",           "atmospheric O2 mixing ratio"),
            ("ignit",   "",           "ignition probability"),
            ("firef",   "",           "effect of fire on VEG"),        
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
            ("silw_relative", "",     "Silicate weathering normalized to present day"),
            ("carbw",   "molC/yr",    "Carbonate weathering"),
            ("carbw_relative", "",    "Carbonate weathering normalized to present day"),
            ("oxw_fac", "",           "Oxidative weathering oxygen-dependent factor"),
            ("oxidw",   "molC/yr",    "Oxidative organic C weathering"),
            ("oxidw_relative", "",    "Oxidative organic C weathering normalized to present day"),
            ("gypw",    "molS/yr",    "Gypsum weathering"),
            ("pyrw",    "molS/yr",    "Pyrite weathering"),
            ("phosw_s", "molP/yr",    "Phosphorus weathering, silicate-associated"),
            ("phosw_c", "molP/yr",    "Phosphorus weathering, carbonate-associated"),
            ("phosw_o", "molP/yr",    "Phosphorus weathering, organic carbon-associated"),
            ("phosw",   "molP/yr",    "Phosphorus weathering, total"),
            ("pland",   "molP/yr",    "Phosphorus weathering to land"),
            ("psea",    "molP/yr",    "Phosphorus weathering to sea"),
            # burial
            ("locb",   "molC/yr",     "Land organic carbon burial"),
            ("locb_delta",   "per mil","Land organic carbon burial d13C"),
            ("eps_eqbw",    "per mil", "approximate atmosphere-runoff d13C")
        ])

    # isotope Types
    CIsotopeType = rj.pars.CIsotope[]
    SIsotopeType = rj.pars.SIsotope[]

    # Add flux couplers
    fluxAtoLand = PB.Fluxes.FluxContribScalar(
        "fluxAtoLand.flux_", ["CO2::$CIsotopeType", "O2"],
        isotope_data=Dict())

    fluxRtoOcean = PB.Fluxes.FluxContribScalar(
        "fluxRtoOcean.flux_", ["DIC::$CIsotopeType", "TAlk", "Ca", "P", "SO4::$SIsotopeType"],
        isotope_data=Dict())

    fluxLandtoSedCrust = PB.Fluxes.FluxContribScalar(
        "fluxLandtoSedCrust.flux_", ["Ccarb::$CIsotopeType", "Corg::$CIsotopeType", "GYP::$SIsotopeType", "PYR::$SIsotopeType"],
        isotope_data=Dict())    

    PB.add_method_do!(
        rj,
        do_land_bergman2004,
        (   
            PB.VarList_namedtuple_fields(fluxAtoLand),
            PB.VarList_namedtuple_fields(fluxRtoOcean),
            PB.VarList_namedtuple_fields(fluxLandtoSedCrust),
            PB.VarList_namedtuple([var_VEG; vars_dep; vars_prop]),
        ),
        p=(CIsotopeType, SIsotopeType),
    )

    return nothing
end


"Calculate rates"
function do_land_bergman2004(
    m::PB.ReactionMethod,
    pars,
    (fluxAtoLand, fluxRtoOcean, fluxLandtoSedCrust, D), 
    cellrange::PB.AbstractCellRange,
    deltat
)
    (Cisotopetype, Sisotopetype) = m.p

    # land biota
    if pars.f_landbiota[] == "Dynamic"
        copse_landbiota_bergman2004( pars, D )
    elseif pars.f_landbiota[] == "Prescribed"
        # VEG supplied as forcing
    else
        error("unknown f_landbiota='$(pars.f_landbiota[])'")
    end

    #%%%%%%%% calculate weathering
    #%%%% Plant effects on weathering
    pCO2PAL_min = 1e-3 # guard against -ve pCO2PAL when using AD
    pCO2PAL_limit = max(D.pCO2PAL[], pCO2PAL_min)

    D.f_preplant[] = copse_f_T(D.TEMP[]) * pCO2PAL_limit^0.5
    D.f_plant[]    = copse_f_T(D.TEMP[]) * ( 2.0*pCO2PAL_limit / (1.0 + pCO2PAL_limit) )^0.4

    D.g_preplant[] = copse_g_T(D.TEMP[]) * pCO2PAL_limit^0.5
    D.g_plant[]    = copse_g_T(D.TEMP[]) * ( 2.0*pCO2PAL_limit / (1.0 + pCO2PAL_limit) )^0.4

    D.VWmin[] = min(D.VEG[]*D.W[], 1.0)

    D.f_co2[] = D.f_preplant[]*(1.0 - D.VWmin[]) + D.f_plant[]*D.VWmin[]
    D.g_co2[] = D.g_preplant[]*(1.0 - D.VWmin[]) + D.g_plant[]*D.VWmin[]

    D.w_plantenhance[] = (pars.k15_plantenhance[] +
                (1.0 - pars.k15_plantenhance[]) * D.W[] * D.VEG[])

    # silicate and carbonate weathering
    D.silw_relative[] = D.UPLIFT[]*D.w_plantenhance[] * D.f_co2[]
    D.silw[] = pars.k_silw[]*D.silw_relative[]

    carbw_fac = D.UPLIFT[]*D.w_plantenhance[]*D.g_co2[]
    if pars.f_carbwC[] == "Cindep"   # Copse 5_14
        D.carbw_relative[] = carbw_fac
    elseif pars.f_carbwC[] == "Cprop"   # A generalization for varying-size C reservoir
        D.carbw_relative[] = carbw_fac*D.C_norm[]
    else
       error("Unknown f_carbwC ", pars.f_carbwC[])
    end
    D.carbw[] = pars.k14_carbw[]*D.carbw_relative[]
    

    #%%%% Oxidative weathering
    # Functional form of oxidative weathering
    if pars.f_oxwO[] == "PowerO2"   # Copse 5_14 base with f_oxw_a = 0.5
        D.oxw_fac[] = D.pO2PAL[]^pars.f_oxw_a[]
    elseif pars.f_oxwO[] == "SatO2"
        D.oxw_fac[] = D.pO2PAL[] /(D.pO2PAL[] + pars.f_oxw_halfsat[])
    else
       error("Unknown f_foxwO ", pars.f_oxwO[])
    end
    # C oxidative weathering
    D.oxidw_relative[]  = D.UPLIFT[]*D.G_norm[]*D.oxw_fac[]
    D.oxidw[]           = pars.k17_oxidw[] * D.oxidw_relative[] 

    # Sulphur weathering
    if pars.enableS[]
        # Gypsum weathering tied to carbonate weathering
        D.gypw[] = pars.k22_gypw[] * D.GYP_norm[] * carbw_fac
        # Pyrite oxidative weathering with same functional form as carbon
        D.pyrw[] = pars.k21_pyrw[] * D.UPLIFT[]*D.PYR_norm[] * D.oxw_fac[]
        # Isotope fractionation of S
        gypw_isotope = @PB.isotope_totaldelta(Sisotopetype, D.gypw[], D.GYP_delta[])
        pyrw_isotope = @PB.isotope_totaldelta(Sisotopetype, D.pyrw[], D.PYR_delta[])
    else
        D.gypw[] = 0.0; gypw_isotope = @PB.isotope_totaldelta(Sisotopetype, 0.0, 0.0)
        D.pyrw[] = 0.0; pyrw_isotope = @PB.isotope_totaldelta(Sisotopetype, 0.0, 0.0)
    end

    # P weathering and delivery to land and sea
    D.phosw_s[] = pars.k10_phosw[] * (2.0/12.0)*D.silw_relative[]   
    D.phosw_c[] = pars.k10_phosw[] * (5.0/12.0)*D.carbw_relative[]
    D.phosw_o[] = pars.k10_phosw[] * (5.0/12.0)*D.oxidw_relative[]
    D.phosw[]   = D.phosw_s[] + D.phosw_c[] + D.phosw_o[]

    D.pland[]   = pars.k11_landfrac[]*D.VEG[]*D.phosw[]
    pland0      = pars.k11_landfrac[]*pars.k10_phosw[]

    D.psea[]    = D.phosw[] - D.pland[]

    #%%%%%%% Burial
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # Land organic carbon burial
    if pars.f_locb[] == "original"
        D.locb[] = pars.k5_locb[] * (D.pland[]/(pland0+eps())) * D.CPland_relative[] # guard against pland0 == 0
    elseif pars.f_locb[] == "Prescribed"
        D.locb[] = pars.k5_locb[]*D.locbpert[]
    else
        error("unknown f_locb='$(pars.f_locb[])")
    end

    ## C isotopes
    ############################################################################
    # Isotope fractionation of organic carbon
    if Cisotopetype <: PB.AbstractIsotopeScalar
        D.locb_delta[]      = D.CO2_delta[] - D.D_P_CO2_locb[]
        # DIC isotopes - need to take account of atmosphere-water fractionation
        # (not critical to get this right, as there is a 'short circuit' atm <-> river -> ocean <-> atm)
        D.eps_eqbw[] = 10.78-0.114*(D.TEMP[] - PB.Constants.k_CtoK)
    end

    locb_isotope       =  @PB.isotope_totaldelta(Cisotopetype, D.locb[], D.locb_delta[])

    oxidw_isotope       = @PB.isotope_totaldelta(Cisotopetype, D.oxidw[], D.G_delta[])

    carbw_isotope       = @PB.isotope_totaldelta(Cisotopetype, D.carbw[], D.C_delta[])

    DICrunoff           = @PB.isotope_totaldelta(Cisotopetype, D.carbw[] + 2*D.silw[], D.CO2_delta[] + D.eps_eqbw[])

    # fluxes
    #########################################################################

    # Atmospheric fluxes 
    fluxAtoLand.CO2[]   += DICrunoff - oxidw_isotope + locb_isotope
    
    fluxAtoLand.O2[]    += D.oxidw[] + 2*D.pyrw[] - D.locb[]

    # Riverine fluxes
    fluxRtoOcean.DIC[]  += DICrunoff + carbw_isotope   

    fluxRtoOcean.TAlk[] += 2*D.silw[] + 2*D.carbw[]    

    fluxRtoOcean.Ca[]   += D.silw[] + D.carbw[] + D.gypw[]

    fluxRtoOcean.P[]    += D.psea[]
    
    fluxRtoOcean.SO4[]  += gypw_isotope + pyrw_isotope

    # Sedimentary fluxes
    fluxLandtoSedCrust.Ccarb[]  += -carbw_isotope

    fluxLandtoSedCrust.Corg[]   += locb_isotope -oxidw_isotope
                
    fluxLandtoSedCrust.GYP[]    += -gypw_isotope

    fluxLandtoSedCrust.PYR[]    += -pyrw_isotope
                
    return nothing
end



function copse_landbiota_bergman2004(pars, D)
    """COPSE_LANDBIOTA COPSE OCT dynamic land vegetation model
    """

    # effect of temp on VEG
    # SD update to BM Matlab - this fixes absence of T limitation
    # visible as slightly high VEG
    D.V_T[] = 1.0 - (((D.TEMP[] - PB.Constants.k_CtoK-25.0)/25.0)^2)

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



function copse_f_T(Tk)
    """ % COPSE_F_T Temperature dependence of silicate weathering from COPSE (Bergman 2004)
        %
        % Tk   - Temperature (Kelvin)
        % f_T  - weathering rate relative to Tzero = 15 C
        %
        % See also COPSE_equations GEOCARB_equations test_copse_weathering
    """

    Tzero = PB.Constants.k_CtoK+15.0
    copse_f_T = exp(0.09*(Tk-Tzero)) * ((1.0 + 0.038*(Tk - Tzero))^0.65)

    return copse_f_T
end

function copse_g_T(Tk)
    """% COPSE_G_T  Temperature dependence of carbonate weathering from COPSE (Bergman 2004)
        %
        % Tk   - Temperature (Kelvin)
        % g_T  - weathering rate relative to Tzero = 15 C
        % NB: will go negative to T <~ 5C !!
        %
        % See also COPSE_equations GEOCARB_equations test_copse_weathering
    """

    Tzero = PB.Constants.k_CtoK+15.0

    g_T = 1.0 + 0.087*(Tk - Tzero)
    if (g_T < 0.0)
        @warn "copse_g_T: negative g_T $(g_T) for Tk $(Tk) resetting to 0.0"
        g_T = 0.0
    end

    return g_T
end


end # module
