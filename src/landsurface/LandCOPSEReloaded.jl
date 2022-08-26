# -*- coding: utf-8 -*-
module LandCOPSEReloaded


import PALEOboxes as PB
using PALEOboxes.DocStrings

import Infiltrator # Julia debugger

"""
    ReactionLandBiota

COPSE Reloaded(2018) land biota. Calculates `VEG`, representing land plant coverage and biomass.

In a scalar (0D) Domain, this functions as in COPSE to calculate a global mean.

In a spatially-resolved Domain, this calculates local properties using local TEMP.

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionLandBiota{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        # land biota
        PB.ParString("f_npp",       "original",
            allowed_values=["original", "noT", "noO2", "noCO2", "bernerCO2", "newCO2", "constant"], 
            description="functional form for VEG limitation by EVO, temperature, O2, CO2"),
        PB.ParDouble("k16_PANtoO",     3.762,      units="",       
            description="For calculating mixing ratio O2 from normalised O2"),
        PB.ParString("f_ignit",     "original",
            allowed_values=["original", "bookchapter", "nofbk"], 
            description="functional form for fire feedback"),
        PB.ParDouble("k_fire",         100.0,      units="",
            description="fire feedback"),
    )

end

function PB.register_methods!(rj::ReactionLandBiota)
            
    # dependencies
    vars_dep_scalar = PB.VarVector(PB.VarDepScalar,
        [
            # Forcings
            ("global.EVO",     "",      "plant evolution scaling"),           
            ("(global.fireforcing)", "",  "optional fireforcing"),
            #       
            # atm
            ("atm.pO2PAL",  "",           "atmospheric oxygen relative to present-day"),
            ("atm.pCO2PAL", "",           "atmospheric pCO2 relative to present-day"),
            ("atm.pCO2atm", "atm",        "atmospheric pCO2"),
        ]
    )

    vars_dep = PB.VarVector(PB.VarDep,
        [
            # Forcings       
            ("global.TEMP",    "K",      "global or local temperature"),
        ]
    )

    # Properties we calculate in do_react
    vars_prop_scalar = PB.VarVector(PB.VarPropScalar,
        [
            # copse_landbiota
            ("V_co2",   "",           "effect of CO2 on VEG"),
            ("V_o2",    "",           "effect of O2 on VEG"),
            ("mrO2",    "",           "atmospheric O2 mixing ratio"),
            ("ignit",   "",           "ignition probability"),
            ("firef",   "",           "effect of fire on VEG"),
        ]
    )

    vars_prop = PB.VarVector(PB.VarProp,
        [
            # copse_landbiota
            ("V_T",     "",           "effect of temp on VEG"),
            ("V_npp",   "",           "full VEG limitation"),
            ("VEG",     "",           "Mass of terrestrial biosphere"),
        ]
    )


    PB.add_method_do!(
        rj,
        do_land_biota,
        (PB.VarList_namedtuple([vars_dep_scalar; vars_dep; vars_prop_scalar; vars_prop]), ),
    )

    return nothing
end

function do_land_biota(m::PB.ReactionMethod, pars, (D, ),  cellrange::PB.AbstractCellRange, deltat)
    # copse_landbiota_reloaded

    # effect of temp on VEG
    for i in cellrange.indices
        D.V_T[i] = max(1 - (( (D.TEMP[i] - PB.Constants.k_CtoK-25)/25 )^2), 0) # defend against low TEMP values eg transients at startup
    end

    # effect of CO2 on VEG
    function V_co2_orig()
        P_atm = D.pCO2atm[]*1e6 
        P_half = 183.6
        P_min = 10.0
        return (P_atm - P_min) / (P_half + P_atm - P_min)
    end    

    # effect of O2 on VEG 
    D.V_o2[] = max(1.5 - 0.5*D.pO2PAL[], 0.0)

    # full VEG limitation
    if pars.f_npp[] == "original"
        D.V_co2[] = V_co2_orig()
        foreach(i -> D.V_npp[i] = 2*D.EVO[]*D.V_T[i]*D.V_o2[]*D.V_co2[], cellrange.indices)
    elseif pars.f_npp[] == "noT"
        D.V_co2[] = V_co2_orig()
        foreach(i -> D.V_npp[i] = 2*D.EVO[]*0.84*D.V_o2[]*D.V_co2[], cellrange.indices)
    elseif pars.f_npp[] == "noO2"
        D.V_co2[] = V_co2_orig()
        D.V_o2[] = NaN
        foreach(i -> D.V_npp[i] = 2*D.EVO[]*D.V_T[i]*1.0*D.V_co2[], cellrange.indices)
    elseif pars.f_npp[] == "noCO2"
        D.V_co2[] = NaN
        foreach(i -> D.V_npp[i] = 2*D.EVO[]*D.V_T[i]*D.V_o2[]*0.5952381, cellrange.indices)
    elseif pars.f_npp[] == "bernerCO2"
        D.V_co2[] = 2*D.pCO2PAL[] / (1.0 + D.pCO2PAL[]) # Berner choice
        foreach(i -> D.V_npp[i] = D.EVO[]*(D.V_T[i]/0.84)*D.V_o2[]*D.V_co2[], cellrange.indices)
    elseif pars.f_npp[] == "newCO2"
        D.V_co2[] = D.pCO2PAL[]^0.5
        foreach(i -> D.V_npp[i] = D.EVO[]*(D.V_T[i]/0.84)*D.V_o2[]*D.V_co2[], cellrange.indices)
    elseif pars.f_npp[] == "constant"
        D.V_co2[] = NaN
        D.V_o2[] = NaN
        foreach(i -> D.V_npp[i] = D.EVO[], cellrange.indices)
    else
        error("Unknown f_npp ", pars.f_npp[])
    end
  

    # fire feedback
    # calculate atmospheric mixing ratio of O2 (for constant atmospheric N etc!)
    #(only used for fire ignition probability)
    D.mrO2[] = D.pO2PAL[]  /   ( D.pO2PAL[]  +pars.k16_PANtoO[] )
    if pars.f_ignit[] == "original"
        D.ignit[] = max(586.2*(D.mrO2[])-122.102 , 0  ) 
    elseif pars.f_ignit[] == "bookchapter"
        # TML new fn. based on data for fuel of 10% moisture content
        D.ignit[] = min(max(48.0*(D.mrO2[])-9.08 , 0  ) , 5  )
    elseif pars.f_ignit[] == "nofbk"
        D.ignit[] = 1.0
    else
        error("Unknown f_ignit ", pars.f_ignit[])
    end

    if isnothing(D.fireforcing)
        fireforcing = 1.0
    else
        fireforcing = D.fireforcing[]
    end
    D.firef[] = fireforcing*pars.k_fire[]/(fireforcing*pars.k_fire[] - 1 + D.ignit[])

    # Mass of terrestrial biosphere
    for i in cellrange.indices
        D.VEG[i] = D.V_npp[i] * D.firef[]
    end

    return nothing
end


"""
    ReactionLandWeatheringRates

COPSE Reloaded(2018) weathering rates.
Calculates weathering rates as a function of temperature, pCO2, pO2 and biosphere `VEG` and biological weathering enhancement `W`.

In a scalar (0D) Domain, this functions as in COPSE to calculate global mean rates.

In a spatially-resolved Domain, this calculates local rates (use 'f_runoff' parameter to choose
an appropriate hydrological parameterisation based on a supplied spatially-varying `runoff (kg m-2 s-1)`).

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionLandWeatheringRates{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        # weathering rates
        PB.ParString("f_act_energies",     "single",
            allowed_values=["single", "split"], 
            description="weathering kinetics for granite and basalt"),
        PB.ParString("f_runoff",     "original",
            allowed_values=["original", "geoclim_2004", "linear_runoff"], 
            description="runoff dependency of weathering"),
        PB.ParDouble("k_runoff_mean", 30e-6, units="kg m-2 s-1",
            description="global mean runoff (f_runoff != original only)"),
        PB.ParDouble("k_runoff_freeze", 0.0, units="K",
            description="local temperature below which weathering is zero (f_runoff != original only)"),
        PB.ParString("f_co2fert",           "original",
            allowed_values=["original", "geocarb3", "off"], 
            description="CO2 and O2 dependence of plant effects on CO2 weathering"),
        PB.ParString("f_oxwO",           "PowerO2",
            allowed_values=["PowerO2", "IndepO2"], 
            description="functional form for oxygen sensitivity of oxidative weathering"),
        PB.ParDouble("f_oxw_a", 0.5,
            description="power law dependence of oxidative weathering oxygen sensitivity"),
        PB.ParDouble("k15_plantenhance",0.15,      units="",
            description="weathering enhancement factor prior to vascular plant colonisation"),
        PB.ParString("f_vegweath",           "original",
            allowed_values=["original", "new", "new2", "newnpp"], 
            description="vegetation dependence of silicate,carbonate weathering"),
    )
end

function PB.register_methods!(rj::ReactionLandWeatheringRates)
        
    # dependencies
    vars_dep_scalar = PB.VarVector(PB.VarDepScalar,
        [
            # Forcings      
            ("global.W",     "",      "biological enhancement of weathering"),   
            #       
            # atm
            ("atm.pO2PAL",  "",           "atmospheric oxygen relative to present-day"),
            ("atm.pCO2PAL", "",           "atmospheric pCO2 relative to present-day"),                
        ]
    )

    vars_dep = PB.VarVector(PB.VarDep,
        [
            # Forcings      
            ("global.TEMP",    "K",      "global or local temperature"),
            ("(runoff)",    "kg m-2 s-1", "runoff"),
            # copse_landbiota     
            ("(V_npp)",   "",           "full VEG limitation"),               
            ("VEG",     "",           "Mass of terrestrial biosphere"), 
        ]
    )

    # Properties we calculate in do_react
    vars_prop_scalar = PB.VarVector(PB.VarPropScalar,
        [
            ("f_CO2_abiotic",   "",         "abiotic CO2 dependence of carbonate silicate weathering"),
            ("f_CO2_biotic",   "",         "biotic CO2 dependence of carbonate silicate weathering"),
            ("oxw_facO",        "",         "oxidative weathering oxygen dependence factor"),
            ("oxw_facOmax",        "",      "maximum possible oxidative weathering factor at high pO2"),            
        ]
    )

    vars_prop = PB.VarVector(PB.VarProp,
        [
            ("f_T_gran",     "",           "granite T dependence from activation energy"),
            ("f_T_bas",     "",           "basalt T dependence from activation energy"),
            ("f_T_ap",     "",           "apatite T dependence from activation energy"),
            ("f_T_runoff",     "",           "runoff silicate temperature dependence"),
            ("g_T_runoff",     "",           "runoff carbonate temperature dependence"),
            ("f_gran",          "",         "granite weathering rate"),
            ("f_bas",          "",         "basalt weathering rate"),
            ("f_ap",          "",         "apatite weathering rate"),
            ("f_carb",          "",         "carbonate weathering rate"), 
            ("VWmin",          "",         "VEG (biosphere) and W (evolutionary) weathering factor"),          
        ]
    )
   
    PB.add_method_do!(
        rj,
        do_land_weathering_rates,
        (PB.VarList_namedtuple([vars_dep_scalar; vars_dep; vars_prop_scalar; vars_prop]), ),       
    )

    return nothing
end


function do_land_weathering_rates(m::PB.ReactionMethod, pars, (D, ),  cellrange::PB.AbstractCellRange, deltat)

    for i in cellrange.indices
        # weathering kinetics 
        if pars.f_act_energies[] == "single"
            # use same activation energy for granite and basalt  
            D.f_T_gran[i] = exp(0.09*(D.TEMP[i]-288.15))
            D.f_T_bas[i] = exp(0.09*(D.TEMP[i]-288.15))
        elseif pars.f_act_energies[] == "split"
            # use different activation energy for granite and basalt  
            D.f_T_gran[i] = exp(0.0724*(D.TEMP[i]-288.15)) # 50 kJ/mol 
            D.f_T_bas[i] = exp(0.0608*(D.TEMP[i]-288.15)) # 42 kJ/mol 
        else
            error("unrecognized f_act_energies ", pars.f_act_energies[])
        end 


        # activation energy for apatite (decide whether to use it in _fluxes) 
        D.f_T_ap[i] = exp(0.0507*(D.TEMP[i]-288.15)) # 35 kJ/mol 

    
        # runoff temperature dependence
        if pars.f_runoff[] == "original"
            D.f_T_runoff[i] = ( max(1 + 0.038*(D.TEMP[i] - 288.15), 0)^0.65 ) # defend against low TEMP values eg transients at startup
            D.g_T_runoff[i] = max(1 + 0.087*(D.TEMP[i] - 288.15),0) # defend against low TEMP values eg transients at startup
        elseif pars.f_runoff[] == "geoclim_2004"
            D.f_T_runoff[i] = (D.TEMP[i] >= pars.k_runoff_freeze[])*max(D.runoff[i], 0.0)/pars.k_runoff_mean[]
            D.g_T_runoff[i] = (D.TEMP[i] >= pars.k_runoff_freeze[])*sqrt(max(D.runoff[i], 0.0)/pars.k_runoff_mean[])
        elseif pars.f_runoff[] == "linear_runoff"
            D.f_T_runoff[i] = (D.TEMP[i] >= pars.k_runoff_freeze[])*max(D.runoff[i], 0.0)/pars.k_runoff_mean[]
            D.g_T_runoff[i] = D.f_T_runoff[i] 
        else
            error("Unknown f_runoff ", pars.f_runoff[])
        end 
    end
 
    # CO2 dependence 
    D.f_CO2_abiotic[] = max(D.pCO2PAL[],0)^0.5
    if pars.f_co2fert[] == "original"
         D.f_CO2_biotic[] = ( 2*max(D.pCO2PAL[],0) / (1 + max(D.pCO2PAL[],0)) )^0.4 
    elseif pars.f_co2fert[] == "geocarb3"
         D.f_CO2_biotic[] = ( 2*max(D.pCO2PAL[],0) / (1 + max(D.pCO2PAL[],0)) )
    elseif pars.f_co2fert[] == "off"
         D.f_CO2_biotic[] = 1.0
    else
        error("Unknown f_co2fert ", pars.f_co2fert[])
    end 
 
 
    # O2 dependence 
    if pars.f_oxwO[] == "PowerO2"   # Copse 5_14 base with f_oxw_a = 0.5
        D.oxw_facO[] = max(D.pO2PAL[], 0)^pars.f_oxw_a[]
        D.oxw_facOmax[] = NaN;
    elseif pars.f_oxwO[] == "SatO2"
        D.oxw_facO[] = max(D.pO2PAL[], 0)/(max(D.pO2PAL[], 0) + pars.f_oxw_halfsat[])
        D.oxw_facOmax[] = 1.0
    elseif pars.f_oxwO[] == "IndepO2"
        D.oxw_facO[] = 1
        D.oxw_facOmax[] = NaN
    # elseif pars.f_oxwO[] == "funcO2"
        # [oxw_O, oxwmax] = pars.f_oxwOfunc( max(pO2PAL, 0) );
        # D.oxw_facO = oxw_O/pars.f_oxwOfunc( 1.0 );
        # D.oxw_facOmax = oxwmax/pars.f_oxwOfunc( 1.0 );
    else
        error("Unknown f_foxwO ", pars.f_oxwO[])
    end
 
 
    # vegetation dependence
    for i in cellrange.indices
        if pars.f_vegweath[] == "original"
            D.VWmin[i] = min(D.VEG[i]*D.W[],1); 
            f_co2_gran = D.f_T_gran[i]*D.f_T_runoff[i]*(D.f_CO2_abiotic[]*(1 - D.VWmin[i]) + D.f_CO2_biotic[]*D.VWmin[i])
            f_co2_bas = D.f_T_bas[i]*D.f_T_runoff[i]*(D.f_CO2_abiotic[]*(1 - D.VWmin[i]) + D.f_CO2_biotic[]*D.VWmin[i])
            f_co2_ap = D.f_T_ap[i]*D.f_T_runoff[i]*(D.f_CO2_abiotic[]*(1 - D.VWmin[i]) + D.f_CO2_biotic[]*D.VWmin[i]) 
            g_co2 = D.g_T_runoff[i]*(D.f_CO2_abiotic[]*(1 - D.VWmin[i]) + D.f_CO2_biotic[]*D.VWmin[i])
            W_plantenhance = (  pars.k15_plantenhance[]  + (1-pars.k15_plantenhance[])*D.W[]*D.VEG[i]  )
            # combining plant and other effects on rates 
            D.f_gran[i] = W_plantenhance * f_co2_gran
            D.f_bas[i] = W_plantenhance *f_co2_bas
            D.f_ap[i] = W_plantenhance * f_co2_ap
            D.f_carb[i] = W_plantenhance* g_co2
        elseif pars.f_vegweath[] == "new"
            D.VWmin[i] = min(D.VEG[i],1);      
            D.f_gran[i] = D.f_T_gran[i]*D.f_T_runoff[i]*((1-D.VWmin[i])*pars.k15_plantenhance[]*D.f_CO2_abiotic[] + D.VEG[i]) ; 
            D.f_bas[i] = D.f_T_bas[i]*D.f_T_runoff[i]*((1-D.VWmin[i])*pars.k15_plantenhance[]*D.f_CO2_abiotic[] + D.VEG[i]) ; 
            D.f_ap[i] = D.f_T_ap[i]*D.f_T_runoff[i]*((1-D.VWmin[i])*pars.k15_plantenhance[]*D.f_CO2_abiotic[] + D.VEG[i]) ; 
            D.f_carb[i] = D.g_T_runoff[i]*((1-D.VWmin[i])*pars.k15_plantenhance[]*D.f_CO2_abiotic[] + D.VEG[i]) ; 
        elseif pars.f_vegweath[] == "new2" 
            D.VWmin[i] = min(D.VEG[i]*D.W[],1); 
            D.f_gran[i] = D.f_T_gran[i]*D.f_T_runoff[i]*((1-D.VWmin[i])*pars.k15_plantenhance[]*D.f_CO2_abiotic[] + D.VEG[i]*D.W[]) ; 
            D.f_bas[i] = D.f_T_bas[i]*D.f_T_runoff[i]*((1-D.VWmin[i])*pars.k15_plantenhance[]*D.f_CO2_abiotic[] + D.VEG[i]*D.W[]) ; 
            D.f_ap[i] = D.f_T_ap[i]*D.f_T_runoff[i]*((1-D.VWmin[i])*pars.k15_plantenhance[]*D.f_CO2_abiotic[] + D.VEG[i]*D.W[]) ; 
            D.f_carb[i] = D.g_T_runoff[i]*((1-D.VWmin[i])*pars.k15_plantenhance[]*D.f_CO2_abiotic[] + D.VEG[i]*D.W[]) ; 
        elseif pars.f_vegweath[] == "newnpp"
            D.VWmin[i] = min(D.V_npp[i]*D.W[],1); 
            D.f_gran[i] = D.f_T_gran[i]*D.f_T_runoff[u]*((1-D.VWmin[i])*pars.k15_plantenhance[]*D.f_CO2_abiotic[] + D.V_npp[i]*D.W[]) ; 
            D.f_bas[i] = D.f_T_bas[i]*D.f_T_runoff[i]*((1-D.VWmin[i])*pars.k15_plantenhance[]*D.f_CO2_abiotic[] + D.V_npp[i]*D.W[]) ; 
            D.f_ap[i] = D.f_T_ap[i]*D.f_T_runoff[i]*((1-D.VWmin[i])*pars.k15_plantenhance[]*D.f_CO2_abiotic[] + D.V_npp[i]*D.W[]) ; 
            D.f_carb[i] = D.g_T_runoff[i]*((1-D.VWmin[i])*pars.k15_plantenhance[]*D.f_CO2_abiotic[] + D.V_npp[i]*D.W[]) ; 
        else
            error("Unknown f_vegweath ", pars.f_vegweath[])
        end
    end

    return nothing
end


"""
    ReactionLandWeatheringFluxes

COPSE Reloaded(2018) weathering fluxes.
Calculates and applies global weathering fluxes and land burial given weathering rates and land areas.

Fluxes are added to flux couplers:
- `fluxAtoLand`:  CO2 and O2 exchange with atmosphere
- `fluxRtoOcean`: riverine fluxes
- `fluxLandtoSedCrust`: sedimentary reservoir weathering and land organic carbon burial

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionLandWeatheringFluxes{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        # silicate weathering
        PB.ParString("f_gran_link_u",       "original",
            allowed_values=["original", "weak", "none"],
            description="granite weathering uplift dependence"),
        PB.ParString("f_bas_link_u",       "yes",
            allowed_values=["yes", "weak", "no"],
            description="basalt weathering uplift dependence"),
        
        PB.ParDouble("k_basfrac",    0.35,   units="",
            description="fraction of silicate weathering that is basaltic at present"),

        PB.ParDouble("k_granw",      NaN,   units="mol/yr",
            description="granite weathering rate"),
        PB.ParDouble("k_basw",      NaN,   units="mol/yr",
            description="basalt weathering rate"),

        # P weathering
        PB.ParString("f_p_kinetics",       "no",
            allowed_values=["no", "yes"],
            description="apatite kinetics independent of host rock"),
        PB.ParString("f_p_apportion",       "no",
            allowed_values=["no", "yes"],
            description="P content vary between granite and basalt"),
        PB.ParDouble("k_P",        2.15,
            description="enrichment P:(Ca+Mg) in granite vs basalt "), 
        PB.ParDouble("k10_phosw",      NaN,        units="mol/yr",
            description="phosphorus weathering"),
        PB.ParDouble("k_Psilw",     2/12,
            description="fraction of P weathering from silicates"),
        PB.ParDouble("k_Pcarbw",     5/12,
            description="fraction of P weathering associated with carbonates"),
        PB.ParDouble("k_Poxidw",     5/12,
            description="fraction of P weathering from organic matter"),
        PB.ParDouble("k_Psedw",      0.0,
            description="fraction of P weathering from other sedimentary rocks"),

        # Carbonate weathering
        PB.ParString("f_carb_link_u",       "yes",
            allowed_values=["yes", "no"],
            description="carbonate weathering uplift dependence"),
        PB.ParString("f_carbwC",       "Cindep",
            allowed_values=["Cindep", "Cprop"],
            description="C (carbonate reservoir) dependence of carbonate weathering"),
        PB.ParDouble("k14_carbw",      13.35e12,   units="mol/yr",
            description="carbonate weathering"),

        # oxidative weathering
        PB.ParString("f_oxwG",       "Gprop",
            allowed_values=["Gindep", "Gprop", "forced"],
            description="G (organic carbon reservoir) dependence of oxidative weathering"),
        PB.ParDouble("k17_oxidw",      NaN,        units="mol/yr",
            description="oxidative org carbon weathering"),

        # Sulphur weathering
        PB.ParBool("enableS",          true,
            description="enable S weathering"),
        PB.ParString("f_pyrweather",    "copse_O2",
            allowed_values=["copse_O2", "copse_noO2", "forced"],
            description="functional form of dependence of pyrite weathering"),
        PB.ParDouble("k21_pyrw",       0.53e12,    units="mol S/yr",
            description="pyrite weathering"),
        PB.ParString("f_gypweather",    "original",
            allowed_values=["original", "alternative", "forced"],
            description="functional form of dependence of gypsum weathering"),
        PB.ParDouble("k22_gypw",       1e12,       units="mol S/yr",
            description="gypsum weathering"),
        PB.ParBool("SRedoxAlk",     false,
            description="true to include -TAlk from pyrite weathering"),

        # burial fluxes
        PB.ParString("f_locb",    "original",
            allowed_values=["original", "Uforced", "coal", "split", "Prescribed"],
            description="functional form of dependence of land organic carbon burial"),
        PB.ParDouble("k5_locb",        4.5e12,     units="mol/yr",
            description="land organic carbon burial (f_locb==Prescribed only)"),
        PB.ParDouble("k11_landfrac",   0.10345,    units="",
            description="fract of weath P buried on land"),
        PB.ParDouble("k_aq",           0.8,       units="",
            description="fraction of locb assumed to occur in aquatic settings (not coals)"),
        PB.ParDouble("CPland0",        1000.0,     units="",
            description="present-day C/P land burial"),

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


function PB.register_methods!(rj::ReactionLandWeatheringFluxes)
 
    # dependencies
    vars_dep = PB.VarVector(PB.VarDepScalar,
        [
            # Forcings
            ("global.tforce",  "yr"   , "time for external forcings"),
            ("global.UPLIFT",  "",      "uplift scaling"),
            ("global.PG",       "",      "paleogeographic forcing"),        
            ("global.RHOSIL",     "",  "silicate weathering enhancement"),
            ("global.RHO",     "",  "carbonate weathering enhancement"),
            ("global.F_EPSILON",     "",  "phosphorus weathering enhancement"),
            ("global.CPland_relative","","land CP burial ratio scaling"),
            ("(global.COAL)",     "",    "coal forcing"),
            # Land areas
            ("GRAN_AREA",     "",   "granite area"),
            ("BA_AREA",     "",     "basalt area"),
            ("CARB_AREA",     "",   "carbonate area"),
            ("(ORG_AREA)",     "",    "organic carbon area"),       
            # Isotopes
            ("(D_P_CO2_locb)", "per mil",  "d13C fractionation between terrestrial organic burial and atmospheric CO2"),
            ("(D_eqbw_CO2)",    "per mil", "d13C fractionation between atmospheric CO2 and fresh water"),
            ("(atm.CO2_delta)", "per mil",  "atmospheric pCO2 delta 13C"),
            # biota
            ("VEG",     "",           "Mass of terrestrial biosphere"),  
            # rates
            ("f_gran",          "",         "granite weathering rate"),
            ("f_bas",          "",         "basalt weathering rate"),
            ("f_ap",          "",         "apatite weathering rate"),
            ("f_carb",          "",         "carbonate weathering rate"),
            ("oxw_facO",        "",         "oxidative weathering oxygen dependence factor"),
            ("oxw_facOmax",        "",      "maximum possible oxidative weathering factor at high pO2"),
            # sedcrust
            ("sedcrust.C_norm",        "",    "Sedimentary carbonate normalized to present day"),
            ("(sedcrust.C_delta)",        "",    "Sedimentary carbonate d13C"),
            ("sedcrust.G_norm",        "",    "Sedimentary organic carbon normalized to present day"),
            ("(sedcrust.G_delta)",        "",    "Sedimentary organic carbon d13C"),
            ("(sedcrust.GYP_norm)",      "",    "Sedimentary gypsum normalized to present day"),
            ("(sedcrust.GYP_delta)",      "",    "Sedimentary gypsum d34S"),
            ("(sedcrust.PYR_norm)",      "",    "Sedimentary pyrite normalized to present day"),
            ("(sedcrust.PYR_delta)",      "",    "Sedimentary pyrite d34S")
        ]
    )

    # Properties we calculate in do_react
    vars_prop = PB.VarVector(PB.VarPropScalar,
        [
            ("granw_relative",    "",    "granite weathering normalized to present"),
            ("granw",    "molC/yr",    "granite weathering"),
            ("basw_relative",    "",    "basalt weathering normalized to present"),
            ("basw",    "molC/yr",    "basalt weathering"),
            ("silw",    "molC/yr",    "total silicate weathering"),
            ("silw_relative", "",     "total silicate weathering normalized to present"),

            ("carbw_relative", "",    "Carbonate weathering normalized to present"),
            ("carbw_fac", "",    "Carbonate weathering normalized to present, without C reservoir dependence"),
            ("carbw",   "molC/yr",    "Carbonate weathering"),
            ("silwcarbw_relative",    "",   "total silicate and carbonate weathering, normalized to present"),
        
            ("orgcw",   "molC/yr",    "organic C erosion"),       
            ("oxidw",   "molC/yr",    "Oxidative organic C weathering"),
            ("gypw",    "molS/yr",    "Gypsum weathering"),
            ("pyrw",    "molS/yr",    "Pyrite weathering"),
            ("granw_ap",    "molC/yr",  "apatite effective granite weathering flux"),      
            ("basw_ap",    "molC/yr",   "apatite effective basalt weathering flux"),
            ("silw_ap",    "molC/yr",   "apatite effective total silicate weathering flux"),
            ("sedw_relative", "",     "Other sediment weathering (sandstones etc) normalized to present"),
            ("phosw_x", "molP/yr",    "Phosphorus weathering, other sediment associated"),
            ("phosw_s", "molP/yr",    "Phosphorus weathering, silicate-associated"),
            ("phosw_c", "molP/yr",    "Phosphorus weathering, carbonate-associated"),
            ("phosw_o", "molP/yr",    "Phosphorus weathering, organic carbon-associated"),
            ("phosw",   "molP/yr",    "Phosphorus weathering, total"),
            ("pland",   "molP/yr",    "Phosphorus weathering to land"),
            ("psea",    "molP/yr",    "Phosphorus weathering to sea"),
            # burial
            ("locb",   "molC/yr",     "Land organic carbon burial"),
            ("locb_delta",   "per mil","Land organic carbon burial d13C")
        ]
    )

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
        do_land_weathering_fluxes,
        (   
            PB.VarList_namedtuple_fields(fluxAtoLand),
            PB.VarList_namedtuple_fields(fluxRtoOcean),
            PB.VarList_namedtuple_fields(fluxLandtoSedCrust),
            PB.VarList_namedtuple([vars_dep; vars_prop]),
        ),
        p=(CIsotopeType, SIsotopeType),
    )

    return nothing
end

function do_land_weathering_fluxes(
    m::PB.ReactionMethod,
    pars,
    (fluxAtoLand, fluxRtoOcean, fluxLandtoSedCrust, D), 
    cellrange::PB.AbstractCellRange,
    deltat
)
    (CIsotopeType, SIsotopeType) = m.p

    # granite weathering  
    if pars.f_gran_link_u[] == "original"
        granw_fac_u = D.UPLIFT[]
    elseif pars.f_gran_link_u[] == "weak"
        granw_fac_u = (D.UPLIFT[]^0.33)
    elseif pars.f_gran_link_u[] == "none"
        granw_fac_u = 1.0
    else
        error("unrecognized f_gran_link_u ", pars.f_gran_link_u[])
    end
    D.granw_relative[]  = granw_fac_u * D.PG[]*D.RHOSIL[] * D.GRAN_AREA[] * D.f_gran[]
    D.granw[]           = D.granw_relative[] * pars.k_granw[]

    # basalt weathering
    if pars.f_bas_link_u[] == "no"
        basw_fac_u = 1.0
    elseif pars.f_bas_link_u[] == "yes"
        basw_fac_u = D.UPLIFT[]
    elseif pars.f_bas_link_u[] == "weak"
        basw_fac_u = D.UPLIFT[]^0.33
    else
        error("unrecognized f_bas_link_u ", pars.f_bas_link_u[])
    end
    D.basw_relative[]   = basw_fac_u * D.PG[]*D.RHOSIL[]*D.BA_AREA[]*D.f_bas[]
    D.basw[]            = D.basw_relative[] * pars.k_basw[]    

    D.silw[] = D.granw[] + D.basw[]
    k_silw = pars.k_granw[] + pars.k_basw[]
    D.silw_relative[] = D.silw[] / k_silw  # define relative rate for use by other modules

    # apatite associated with silicate weathering
    # does it follow its own kinetics? (or that of the host rock)
    if pars.f_p_kinetics[] == "no"
        D.granw_ap[] = D.granw[]
        D.basw_ap[] = D.basw[]
    elseif pars.f_p_kinetics[] == "yes"
        D.granw_ap[] = D.granw[] * (D.f_ap[] / D.f_gran[])
        D.basw_ap[] = D.basw[] * (D.f_ap[] / D.f_bas[])
    else
        error("unrecognized f_p_kinetics ", pars.f_p_kinetics[])
    end
    
    # is the P content assumed to vary between granite and basalt?
    if pars.f_p_apportion[] == "no"
        D.silw_ap[] = D.granw_ap[] + D.basw_ap[]
    elseif pars.f_p_apportion[] == "yes"
        D.silw_ap[] = (pars.k_P[]*D.granw_ap[] + D.basw_ap[])/(pars.k_P[]*(1-pars.k_basfrac[])+pars.k_basfrac[])
    else
        error("unrecognized f_p_apportion ", pars.f_p_apportion[])
    end
    
    D.phosw_s[] = D.F_EPSILON[]*pars.k10_phosw[]*pars.k_Psilw[]*(D.silw_ap[]/(k_silw + eps()))  # trap 0/0 if k_silw = 0
    

    # carbonate weathering 
    if pars.f_carb_link_u[] == "yes"
        carbw_fac_u = D.UPLIFT[]
    elseif pars.f_carb_link_u[] == "no"
        carbw_fac_u = 1.0
    else
        error("Unknown f_carb_link_u ", pars.f_carb_link_u[])
    end
    D.carbw_fac[] = carbw_fac_u * D.PG[]*D.RHO[]* D.CARB_AREA[] * D.f_carb[]
  
    if pars.f_carbwC[] == "Cindep"   # Copse 5_14
        carbw_fac_C = 1.0
    elseif pars.f_carbwC[] == "Cprop"    # A generalization for varying-size C reservoir
        carbw_fac_C = D.C_norm[]
    else
        error("Unknown f_carbw ", pars.f_carbwC[])
    end
    D.carbw_relative[] = carbw_fac_C * D.carbw_fac[]
    D.carbw[] = D.carbw_relative[] * pars.k14_carbw[] 
    
    # Define a normalized weathering rate for use by tracers etc.
    D.silwcarbw_relative[] = ( ( D.silw[] + D.carbw[] ) / ( k_silw + pars.k14_carbw[] ) )

    # Oxidative weathering

    # C oxidative weathering
    # not affected by Dforce.PG in G3 paper...
    if pars.f_oxwG[] == "Gindep"
        oxw_fac_G = 1.0
    elseif pars.f_oxwG[] == "Gprop"
        oxw_fac_G = D.G_norm[]    
    elseif pars.f_oxwG[] == "forced"
        oxw_fac_G = D.G_norm[]  * D.ORG_AREA[]
    else
        error("Unknown f_oxwG ", pars.f_oxwG[])
    end
    D.orgcw[] = oxw_fac_G * pars.k17_oxidw[]*D.UPLIFT[] * D.oxw_facOmax[]
    D.oxidw[] = oxw_fac_G * pars.k17_oxidw[]*D.UPLIFT[] * D.oxw_facO[]

    # Sulphur weathering
    if pars.enableS[]
        # Gypsum weathering
        # not tied to Dforce.PG in G3 paper...
        if pars.f_gypweather[] == "original" # Gypsum weathering tied to carbonate weathering
            D.gypw[] = pars.k22_gypw[] * D.GYP_norm[]*D.carbw_fac[]
        elseif pars.f_gypweather[] == "alternative" # independent of carbonate area
            D.gypw[] = pars.k22_gypw[]*D.GYP_norm[]*D.UPLIFT[]*D.PG[]*D.RHO[]*D.f_carb[]
        elseif pars.f_gypweather[] == "forced" # dependent on evaporite area
            D.gypw[] = pars.k22_gypw[]*D.GYP_norm[]*D.EVAP_AREA[]*D.UPLIFT[]*D.PG[]*D.RHO[]*D.f_carb[]
        else
            error("unknown f_gypweather ", pars.f_gypweather[])
        end

        # Pyrite oxidative weathering 
        # not tied to Dforce.PG in G3 paper...
        if pars.f_pyrweather[] == "copse_O2"
            # with same functional form as carbon
            D.pyrw[] = pars.k21_pyrw[]*D.UPLIFT[]*D.PYR_norm[]*D.oxw_facO[]
        elseif pars.f_pyrweather[] == "copse_noO2"    # independent of O2
            D.pyrw[] = pars.k21_pyrw[]*D.UPLIFT[]*D.PYR_norm[]
        elseif pars.f_pyrweather[] == "forced" # forced by exposed shale area
            D.pyrw[] = pars.k21_pyrw[]*D.UPLIFT[]*D.SHALE_AREA[]*D.PYR_norm[]
        else
            error("unknown f_pyrweather ", pars.f_pyrweather[])
        end
        # Isotope fractionation of S
        gypw_isotope = @PB.isotope_totaldelta(SIsotopeType, D.gypw[], D.GYP_delta[])
        pyrw_isotope = @PB.isotope_totaldelta(SIsotopeType, D.pyrw[], D.PYR_delta[])
    else
        D.gypw[] = 0.0; gypw_isotope = @PB.isotope_totaldelta(SIsotopeType, 0.0, 0.0)
        D.pyrw[] = 0.0; pyrw_isotope = @PB.isotope_totaldelta(SIsotopeType, 0.0, 0.0)
    end

    # P weathering 
    # D.phosw_s is defined above

    # Introduction of P weathering flux from sandstones etc
    D.sedw_relative[] = D.UPLIFT[]*D.PG[]*D.RHO[]*D.VEG[]
    D.phosw_x[] = D.F_EPSILON[]*pars.k10_phosw[]*pars.k_Psedw[]*D.sedw_relative[]
    
    D.phosw_c[] = D.F_EPSILON[]*pars.k10_phosw[]*pars.k_Pcarbw[]*(D.carbw[]/(pars.k14_carbw[] + eps()))
    D.phosw_o[] = D.F_EPSILON[]*pars.k10_phosw[]*pars.k_Poxidw[]*(D.oxidw[]/(pars.k17_oxidw[] + eps()))  # trap 0/0 if k17_oxidw = 0
    D.phosw[]   = D.phosw_s[] + D.phosw_c[] + D.phosw_o[] + D.phosw_x[]


    #%%%%%%% Burial
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # Land organic carbon burial
    if     pars.f_locb[] == "original"
        D.pland[] = pars.k11_landfrac[]*D.VEG[]*D.phosw[]
    elseif pars.f_locb[] == "Uforced"
        # Uplift/erosion control of locb
        D.pland[] = pars.k11_landfrac[]*D.UPLIFT[]*D.VEG[]*D.phosw[]
    elseif pars.f_locb[] == "coal"
        # Coal basin forcing of locb
        D.pland[] = pars.k11_landfrac[]*D.VEG[]*D.phosw[]*(pars.k_aq[]+(1-pars.k_aq[])*D.COAL[])
    elseif pars.f_locb[] == "split"
        # Separating aquatic and coal basin components of locb
        D.pland[] = pars.k11_landfrac[]*D.VEG[]*D.phosw[]*(pars.k_aq[]*D.UPLIFT[]+(1-pars.k_aq[])*D.COAL[])
    elseif pars.f_locb[] == "Prescribed"
        # locb forced
        D.pland[] = pars.k11_landfrac[]*D.phosw[]
    else
        error("unknown f_locb ", pars.f_locb[])
    end

    D.psea[] = D.phosw[] - D.pland[]

    # Land organic carbon burial
    if pars.f_locb[] != "Prescribed"
        D.locb[] = D.pland[]*pars.CPland0[]*D.CPland_relative[]           
    end
   

    ## C isotopes
    ############################################################################
    # Isotope fractionation of organic carbon
    if CIsotopeType <: PB.AbstractIsotopeScalar
        D.locb_delta[]      = D.CO2_delta[] - D.D_P_CO2_locb[]        
    end
    locb_isotope       =  @PB.isotope_totaldelta(CIsotopeType, D.locb[], D.locb_delta[])

    oxidw_isotope       = @PB.isotope_totaldelta(CIsotopeType, D.oxidw[], D.G_delta[])

    carbw_isotope       = @PB.isotope_totaldelta(CIsotopeType, D.carbw[], D.C_delta[])

    # DIC isotopes - need to take account of atmosphere-water fractionation
    # (not critical to get this right, as there is a 'short circuit' atm <-> river -> ocean <-> atm)
    DICrunoff           = @PB.isotope_totaldelta(CIsotopeType, D.carbw[] + 2*D.silw[], D.CO2_delta[] + D.D_eqbw_CO2[])

    # fluxes
    #########################################################################

    # Atmospheric fluxes 
    fluxAtoLand.CO2[]   += DICrunoff - oxidw_isotope + locb_isotope
    
    fluxAtoLand.O2[]    += D.oxidw[] + 2*D.pyrw[] - D.locb[]

    # Riverine fluxes
    fluxRtoOcean.DIC[]  += DICrunoff + carbw_isotope   

    fluxRtoOcean.TAlk[] += 2*D.silw[] + 2*D.carbw[]
    if pars.SRedoxAlk[]
        fluxRtoOcean.TAlk[] += -2*D.pyrw[]
    end

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


"""
    set_steady_state(landfluxes::ReactionLandWeatheringFluxes, mocb, mtotpb, C_degass, S_degass, sfw)
    
Set weathering parameters (silicate, carbon oxidation, phosphorus) for modern steady-state for COPSE Reloaded model configuration.

Arguments:
- `mocb`:       mol C yr-1 marine organic carbon burial
- `mtotp`:      mol P yr-1 marine total phosphorus burial ([] to disable P balance)
- `C_degass`:   mol C yr-1 total C (Corg + Ccarb) degassing    
- `S_degass`:  'mol S yr-1 total S degassing (to include degassing in alk balance)
- `sfw`:        mol C yr-1 seafloor weathering 
"""
function set_steady_state(
    landfluxes::ReactionLandWeatheringFluxes, 
    mocb,
    mtotpb,
    C_degass_carb,
    C_degass_org,
    S_degass,
    sfw
)
    @info "LandCOPSEReloaded.set_steady_state:"

    # oxidw
    PB.setvalue!(landfluxes.pars.k17_oxidw, mocb + landfluxes.pars.k5_locb[] - C_degass_org)
    @info "    set k17_oxidw=$(landfluxes.pars.k17_oxidw[]) (mol yr-1) to balance degassing and carbon burial"

    # silicate weathering to balance degassing +/- organic C cycle
    k_silw  = (landfluxes.pars.k17_oxidw[] - (mocb + landfluxes.pars.k5_locb[])
                + C_degass_carb + C_degass_org - sfw)
    k_silw  += S_degass # SD include sulphur degassing ("S" -> CaSO4 ?) 
    # (TODO confirm correct S to Alk factor: presumably should also be a pyr/gyp weathering-burial balance term cf Corg)
    
    PB.setvalue!(landfluxes.pars.k_granw, k_silw*(1-landfluxes.pars.k_basfrac[]))
    PB.setvalue!(landfluxes.pars.k_basw, k_silw*(landfluxes.pars.k_basfrac[]))
    @info "    set k_granw=$(landfluxes.pars.k_granw[]), k_basw=$(landfluxes.pars.k_basw[]) (mol yr-1) "*
        " to balance degassing and organic carbon fluxes"

    # P weathering to balance P burial
    if !isnothing(mtotpb)
        # run.tm.pars.k10_phosw      = ( run.tm.pars.k2_mocb/run.tm.pars.CPsea0  +run.tm.pars.k7_capb +run.tm.pars.k6_fepb )  / (1-run.tm.pars.k11_landfrac) ;
        # TL alternative approach to avoid redundancy and derive k11
        PB.setvalue!(
            landfluxes.pars.k10_phosw,
            mtotpb  + (landfluxes.pars.k5_locb[]/landfluxes.pars.CPland0[])
        )        
        PB.setvalue!(
            landfluxes.pars.k11_landfrac,
            (landfluxes.pars.k5_locb[]/landfluxes.pars.CPland0[])/landfluxes.pars.k10_phosw[]
        )
        @info "    set k10_phosw=$(landfluxes.pars.k10_phosw[]) (mol P yr-1), "*
            "k11_landfrac=$(landfluxes.pars.k11_landfrac[]) to balance P burial"
    end          
end



end # module
