# -*- coding: utf-8 -*-
module OceanCOPSE


import PALEOboxes as PB
using PALEOboxes.DocStrings

import PALEOcopse

"""
    ReactionOceanCOPSE

COPSE Bergman(2004), COPSE Reloaded Lenton etal (2018) 0D ocean

Fluxes are added to flux couplers:
- `fluxOceanBurial`: ocean burial fluxes

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionOceanCOPSE{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        ##### model constants (Table 4 Bergman etal 2004, consts.hh in COPSE 5_14 C code)
        PB.ParDouble("k1_oxfrac",  0.86,       units="",       description="initial oxic fraction"),
        PB.ParDouble("k2_mocb",    4.5e12,     units="mol/yr", description="ocean organic carbon burial"),
        PB.ParDouble("k3_nfix",    8.75e12,    units="mol/yr", description="nitrogen fixation"),
        PB.ParDouble("k4_denit",   4.3e12,     units="mol/yr", description="denitrification"),
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


function PB.register_methods!(rj::ReactionOceanCOPSE)

  
    # TODO - not the right place for this ? 
    # define a 1 cell grid
    rj.domain.grid = PB.Grids.UnstructuredVectorGrid(ncells=1)
    PB.Grids.set_subdomain!(rj.domain.grid, "oceansurface", PB.Grids.BoundarySubdomain([1]), true)
    PB.Grids.set_subdomain!(rj.domain.grid, "oceanfloor",   PB.Grids.BoundarySubdomain([1]), true)

    # isotope Types
    CIsotopeType = rj.pars.CIsotope.v
    if rj.pars.enableS.v
        SIsotopeType = rj.pars.SIsotope.v
    else
        SIsotopeType = PB.ScalarData
    end

    state_varnames = [
        ("P",                   "mol P",    "Marine phosphorus"),         
        ("O",                   "mol O2",   "Atm-ocean oxygen"),
        ("(DIC::$CIsotopeType)",    "mol C",    "ocean inorganic carbon"),
        ("(CAL)",               "mol Ca",   "Marine calcium"),
    ]
    if rj.pars.enableS.v
        push!(state_varnames, ("S::$SIsotopeType",      "mol S",    "Marine sulphate"))
    end
    if rj.pars.f_ncycle.v
        push!(state_varnames, ("N",               "mol N",    "Marine nitrogen"))
    end
    vars_res, vars_sms, vars_dep_res = PB.Reservoirs.ReservoirLinksVector(
        Dict(), state_varnames
    )

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
             # Special-case a dependency on riverine alk flux for carbonate burial
            ("fluxRtoOcean.flux_TAlk",  "mol yr-1"   , "riverine alkalinity flux"),
            ("(oceanfloor.sfw_relative)",    "",      "seafloor weathering relative to present"),
        ]
    )
            
    # Properties we calculate in do_react
    varnames_prop_react = [
        # copse_marinebiota
        ("Pconc",   "molP/kgsw",  "marine P conc"),
        ("newp",    "",           "marine new production"),
        ("ANOX",    "",           "marine anoxic fraction"),
        
        # C isotopes
        ("mccb_delta", "per mil", "D13C fractionation of marine calcium carbonate burial"),
        ("mocb_delta", "per mil", "D13C fractionation of marine organic carbon burial"),
        ("D_mccb_A",    "per mil", "d13C fractionation between marine carbonate burial and global DIC+CO2"),    

        # burial
        ("mccb",   "molC/yr",     "Carbonate burial"),
        ("CPsea",   "",           "marine C:P ratio"),
        ("mocb",   "molC/yr",     "Marine organic carbon burial"),       
    ]
    if rj.pars.enableS.v
        push!(varnames_prop_react, ("D_mpsb", "per mil",     "D34S fractionation pyrite - marine sulphate"))
    end
    if rj.pars.f_ncycle.v
        push!(varnames_prop_react,
            ("Nconc",   "molN/kgsw",  "marine N conc"),
            ("nfix",    "molN/yr",    "marine nitrogen fixation"),
            ("denit",   "molP/yr",    "marine denitrification"),
            ("monb",   "molN/yr",     "Marine organic N burial"),
        )
    end
    vars_prop_react = PB.VarVector(PB.VarPropScalar, varnames_prop_react)

    # Define flux coupler (a NamedTuple of Variables)
    burial_fluxnames = [   
        "Corg::$CIsotopeType", "Ccarb::$CIsotopeType", 
        "Porg", "Pauth", "PFe", "P", 
    ]
    if rj.pars.enableS.v
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
function do_react(m::PB.ReactionMethod,  (fluxOceanBurial, S, D), cellrange::PB.AbstractCellRange, deltat)
    pars = m.reaction.pars
    (CIsotopeType, SIsotopeType) = m.p

    ###################
    #  Marine biota
    ###################

    copse_marinebiota(pars, D.tforce[], S, D )

    #################
    # Isotopes
    #################
    if CIsotopeType <: PB.AbstractIsotopeScalar
        # delta of marine carbonate burial
        D.mccb_delta[] = D.DIC_delta[] + D.D_mccb_DIC[]
       
        # delta of marine organic carbon burial
        D.mocb_delta[] = D.mccb_delta[] - D.D_B_mccb_mocb[]

        # for diagnostic output only
        D.D_mccb_A[] = D.D_oceanDIC_A[] + D.D_mccb_DIC[]
    end

    if SIsotopeType <: PB.AbstractIsotopeScalar
        # Pyrite sulphur isotope fractionation relative to sulphate and gypsum
        if pars.f_sisotopefrac.v == "fixed"
            D.D_mpsb[] = 35.0
        elseif pars.f_sisotopefrac.v == "copse_O2"
            D.D_mpsb[] = 35.0*D.O_norm[]
        else
            error("unknown f_sisotopefrac ", pars.f_sisotopefrac.v)
        end
    end

    ###########
    # Burial
    ###########

    #%%%%% Reduced C species burial
    # Marine organic carbon burial
    if      pars.f_mocb.v == "original"
        D.mocb[]                = pars.k2_mocb.v * (D.newp[]/pars.newp0.v)^pars.f_mocb_b.v
    elseif  pars.f_mocb.v == "Uforced"
        D.mocb[]                = pars.k2_mocb.v * (D.newp[]/pars.newp0.v)^pars.f_mocb_b.v * D.UPLIFT[]
    elseif  pars.f_mocb.v == "O2dep"
        D.mocb[]                = pars.k2_mocb.v * (D.newp[]/pars.newp0.v)^pars.f_mocb_b.v * 2.1276*exp(-0.755*D.pO2PAL[])
    elseif  pars.f_mocb.v == "both"
        D.mocb[]                = pars.k2_mocb.v * (D.newp[]/pars.newp0.v)^pars.f_mocb_b.v * D.UPLIFT[] * 2.1276*exp(-0.755*D.pO2PAL[])
    else
        error("unknown f_mocb ", pars.f_mocb.v)
    end
    mocb_isotope            = @PB.isotope_totaldelta(CIsotopeType, D.mocb[], D.mocb_delta[])
    
    # update tendencies and external fluxes
    S.O_sms[]               += D.mocb[]

    S.DIC_sms[]             += -mocb_isotope
    fluxOceanBurial.Corg[]  += mocb_isotope    

     #%%%%% Oxidised C species burial    
    D.mccb[]                = 0.5*D.flux_TAlk[]  # alkalinity balance
    mccb_isotope            = @PB.isotope_totaldelta(CIsotopeType, D.mccb[], D.mccb_delta[])
    
    PB.add_if_available(S.CAL_sms, -D.mccb[])
    S.DIC_sms[]             += -mccb_isotope
    fluxOceanBurial.Ccarb[] += mccb_isotope

    #%%%% CP ratio
    if   pars.f_CPsea.v == "Fixed"
        D.CPsea[] = pars.CPsea0.v
    elseif pars.f_CPsea.v == "VCI"  # NB typo in Bergman (2004) has dependency reversed
        D.CPsea[] = (pars.f_CPsea_VCI_oxic.v*pars.f_CPsea_VCI_anoxic.v /
            ((1.0-D.ANOX[])*pars.f_CPsea_VCI_anoxic.v + D.ANOX[]*pars.f_CPsea_VCI_oxic.v))
    else
        error("unrecognized f_CPsea ", pars.f_CPsea.v)
    end

    # Marine organic P burial
    mopb                    = (D.mocb[]/D.CPsea[])
    # Marine carbonate-associated P burial
    if      pars.f_capb.v == "original"
        capb                = pars.k7_capb.v * ((D.newp[]/pars.newp0.v)^pars.f_mocb_b.v)
    elseif  pars.f_capb.v == "redox"
        capb                = (pars.k7_capb.v * ((D.newp[]/pars.newp0.v)^pars.f_mocb_b.v)
                                    *(0.5+0.5*(1.0-D.ANOX[])/pars.k1_oxfrac.v))
    else
        error("unknown f_capb ", pars.f_capb.v)
    end
    # Marine Fe-sorbed P burial NB: COPSE 5_14 uses copse_crash to limit at low P
    if      pars.f_fepb.v == "original"
        fepb = ((pars.k6_fepb.v/pars.k1_oxfrac.v) * (1.0 - D.ANOX[]) *
                PALEOcopse.COPSE.copse_crash(D.P_norm[], "fepb", D.tforce[]) )
    elseif  pars.f_fepb.v == "Dforced"
        fepb = D.DEGASS[]*((pars.k6_fepb.v/pars.k1_oxfrac.v) * (1.0 - D.ANOX[]) *
                PALEOcopse.COPSE.copse_crash(D.P_norm[], "fepb", D.tforce[]) )
    elseif  pars.f_fepb.v == "sfw"
        fepb = D.sfw_relative[]*((pars.k6_fepb.v/pars.k1_oxfrac.v) * (1.0 - D.ANOX[]) *
                                PALEOcopse.COPSE.copse_crash(D.P_norm[], "fepb", D.tforce[]) )
    elseif  pars.f_fepb.v == "pdep"
        fepb = pars.k6_fepb.v/pars.k1_oxfrac.v*(1.0-D.ANOX[])*D.P_norm[]
    else
        error("unknown f_fepb ", pars.f_fepb.v)
    end

    totpb                   = mopb + capb + fepb    

    # update tendencies and external fluxes
    S.P_sms[] += -totpb
    fluxOceanBurial.Porg[]  += mopb
    fluxOceanBurial.Pauth[] += capb
    fluxOceanBurial.PFe[]   +=  fepb
    fluxOceanBurial.P[]     +=  totpb

    if pars.f_ncycle.v
        # Marine organic nitrogen burial
        D.monb[]                = D.mocb[]/pars.CNsea0.v
        # update tendencies, no external flux
        S.N_sms[]               += D.nfix[] - D.denit[] - D.monb[]
    end

    # Marine sulphur burial
    if pars.enableS.v
        # Marine gypsum sulphur burial

        mgsb                    = pars.k_mgsb.v * D.S_norm[] * D.CAL_norm[]
        mgsb_isotope            = @PB.isotope_totaldelta(SIsotopeType, mgsb, D.S_delta[])    

        # Marine pyrite sulphur burial
        if pars.f_pyrburial.v == "copse_noO2"  # dependent on sulphate and marine carbon burial
            mpsb                    = pars.k_mpsb.v*D.S_norm[]*(D.mocb[]/pars.k2_mocb.v)
        elseif pars.f_pyrburial.v == "copse_O2"   # dependent on oxygen, sulphate, and marine carbon burial
            mpsb                    = pars.k_mpsb.v*D.S_norm[]/D.O_norm[]*(D.mocb[]/pars.k2_mocb.v)
        else
            error("unknown f_pyrburial ", pars.f_pyrburial.v)
        end
        mpsb_isotope            = @PB.isotope_totaldelta(SIsotopeType, mpsb, D.S_delta[] - D.D_mpsb[])

        # Update tendencies and external fluxes
        S.S_sms[]               += -(mgsb_isotope + mpsb_isotope)
        
        PB.add_if_available(S.CAL_sms, -mgsb)   # CAL is optional (eg not used in COPSE reloaded configs)
        S.O_sms[]               += 2*mpsb

        fluxOceanBurial.GYP[]   += mgsb_isotope
        fluxOceanBurial.PYR[]   += mpsb_isotope
    end
    
    return nothing
end

function copse_marinebiota(pars, tmodel, S, D )
    """COPSE_MARINEBIOTA COPSE marine ecosystem model
    """

    # convert marine nutrient reservoir moles to micromoles/kg concentration
    D.Pconc[] = D.P_norm[] * 2.2  
    #  clunky way to get normalization values
    # fails as Reaction doesn't see norm_value for variable it is linked to
    # P0 = S.P.norm_value
    # N0 = S.N.norm_value
    P0 = S.P[]/D.P_norm[]

    if pars.f_ncycle.v
        D.Nconc[] = D.N_norm[] * 30.9 ;
        N0 = S.N[]/D.N_norm[]
    end


    if pars.f_ncycle.v
        D.newp[] = 117.0 * min(D.Nconc[]/16.0, D.Pconc[])
    else
        D.newp[] = 117.0 * D.Pconc[]
    end

    #%%%%%% OCEAN ANOXIC FRACTION
    if      pars.f_anoxia.v == "original"
        D.ANOX[] = max( 1.0 - pars.k1_oxfrac.v*D.pO2PAL[] * pars.newp0.v/D.newp[], 0.0 )
    elseif  pars.f_anoxia.v == "newanoxia"
        D.ANOX[] = 1/(1 + exp(-pars.k_logistic.v*(pars.k_uptake.v*(D.newp[]/pars.newp0.v)-D.pO2PAL[]) ) )
    else
        error("unknown f_anoxia ", pars.f_anoxia.v)
    end

    #%%%%% nitrogen cycle
    if pars.f_ncycle.v
        if (S.N[]/16.0) < S.P[] 
            D.nfix[] = pars.k3_nfix.v *( ( S.P[] - S.N[]/16.0 ) / ( P0  - N0/16.0) )^pars.f_nfix_power.v
        else
            if   pars.f_nfix_nreplete.v == "Off" # Surely more defensible ?
                D.nfix[] = 0.0
            elseif pars.f_nfix_nreplete.v == "Sign" # SD - COPSE 5_14 C code has this (?!)
                D.nfix[] = pars.k3_nfix.v *( - ( S.P[] - S.N[]/16.0 ) / ( P0  - N0/16.0)  )^pars.f_nfix_power.v
                print("COPSE_equations -ve nfix (check pars.f_nfix_nreplete) tmodel ", tmodel)
            else
                error("unrecognized f_nfix_nreplete ", pars.f_nfix_nreplete.v)
            end
        end
        # Denitrification 
        if      pars.f_denit.v == "original" # NB: COPSE 5_14 uses copse_crash to limit at low N
            D.denit[] = (pars.k4_denit.v * (1.0 + D.ANOX[] / (1.0 - pars.k1_oxfrac.v) ) *
                        PALEOcopse.COPSE.copse_crash(D.N_norm[], "denit", tmodel))
        elseif  pars.f_denit.v == "new" # introduce dependency on [NO3] throughout
            D.denit[] = pars.k4_denit.v * (1.0 + D.ANOX[] / (1.0 - pars.k1_oxfrac.v)) * D.N_norm[]
        else
            error("unknown f_denit ", pars.f_denit.v)
        end
    end

    return nothing
end



end # module
