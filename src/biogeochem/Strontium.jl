"""
    Strontium

Strontium system with 87Sr/86Sr.

Based on [Francois1992](@cite). This implementation is described in  [Mills2014b](@cite), [Lenton2018](@cite).

A Sr model configuration should contain the following Reservoirs and Reactions:

|Domain name         |Reservoirs     |  Reactions                           |
|:-------------------|:--------------|:-------------------------------------|
|sedcrust            |Sr_sed         |`ReactionSrMantleCrust`, `ReactionSrSed`| 
|land                |               |`ReactionSrLand`  | 
|ocean               |Sr               |                |
|oceanfloor          |               |`ReactionSrOceanfloor`  | 
"""
module Strontium

import PALEOboxes as PB

import Infiltrator # Julia debugger

"""
    ReactionSrMantleCrust

Calculate strontium isotope composition of mantle, old (granite), and new (basalt) igneous rocks as a function of time,
given present day isotopic composition and an initial uniform value at the formation of the Earth.
Assumes all subsequent evolution was due to Rb decay with different Rb/Sr ratios (ie no exchange between these mantle and rock reserovirs).
"""
Base.@kwdef mutable struct ReactionSrMantleCrust{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("Sr_delta_0", 0.69898, units="",
            description="uniform d78Sr isotopic fractionation at formation of Earth"),
        PB.ParDouble("lambda_Rb", 1.4e-11, units="yr-1",
            description="Rb decay rate to 87Sr"),

        PB.ParDouble("Sr_old_ig_delta_present", NaN, units="",
            description="present-day d78Sr isotopic composition of old igneous rocks (granite)"),
        PB.ParDouble("Sr_new_ig_delta_present", NaN, units="",
            description="present-day d78Sr isotopic composition of new igneous rocks (basalt)"),
        PB.ParDouble("Sr_mantle_delta_present", NaN, units="",
            description="present-day d78Sr isotopic composition of mantle"),
    )
end

function PB.register_methods!(rj::ReactionSrMantleCrust)

    vars = [
        PB.VarDepScalar("global.tforce", "yr",  "historical time at which to apply forcings, present = 0 yr"),
        PB.VarPropScalar("Sr_old_ig_delta",  "", "d87Sr old igneous rocks (granite)"),
        PB.VarPropScalar("Sr_new_ig_delta",  "", "d87Sr new igneous rocks (basalt)"),
        PB.VarPropScalar("Sr_mantle_delta",  "", "d87Sr mantle"),
    ]

    PB.add_method_do!(
        rj, 
        do_Sr_mantle_crust,
        (PB.VarList_namedtuple(vars), ),
    )

    return nothing
end

"calculate d87Sr at time tyr relative to present-day (-ve) given present-day d87Sr,
 assuming all evolution is due to Rb decay"
function d87Sr_tyr(rj::ReactionSrMantleCrust, d87Sr_present, tyr)

    # calculate Rb to Sr ratio at present day required to give d87Sr_present
    RbSr = (d87Sr_present - rj.pars.Sr_delta_0.v) / 
                (1.0 - exp(-rj.pars.lambda_Rb.v*PB.Constants.age_present_yr))

    tforwards = PB.Constants.age_present_yr + tyr # get time since formation of Earth at tyr

    d87Sr   = rj.pars.Sr_delta_0.v + RbSr*(1.0 - exp(-rj.pars.lambda_Rb.v*tforwards))

    return d87Sr
end

function do_Sr_mantle_crust(
    m::PB.ReactionMethod,
    (vars, ),
    cellrange::PB.AbstractCellRange,
    deltat
)
    rj = m.reaction

    vars.Sr_old_ig_delta[] = d87Sr_tyr(rj, rj.pars.Sr_old_ig_delta_present.v, vars.tforce[])
    vars.Sr_new_ig_delta[] = d87Sr_tyr(rj, rj.pars.Sr_new_ig_delta_present.v, vars.tforce[])
    vars.Sr_mantle_delta[] = d87Sr_tyr(rj, rj.pars.Sr_mantle_delta_present.v, vars.tforce[])

    return nothing
end


"""
    ReactionSrSed

Calculate contribution to rate of change of a sedimentary `Sr_sed` reservoir due to metamorphic loss, and Rb decay.
"""
Base.@kwdef mutable struct ReactionSrSed{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("lambda_Rb", 1.4e-11, units="yr-1",
            description="Rb decay rate to 87Sr"),
        PB.ParDouble("sediment_RbSr", 0.5,
            description="present-day sediment Rb:Sr ratio"),

        PB.ParString("f_Sr_metam", "alternative", allowed_values=["original", "alternative"],
            description="functional form for sediment Sr metamorphic loss"),
        PB.ParDouble("k_Sr_metam", 13e9, units="mol yr-1",
            description="rate of metamorphic loss of Sr from sedimentary reservoir"),
    )
end

function PB.register_methods!(rj::ReactionSrSed)
    _, SrIsotopeType = PB.split_nameisotope("::SrIsotope", rj.external_parameters)

    vars = [
        PB.VarDepScalar("global.tforce", "yr",  "historical time at which to apply forcings, present = 0 yr"),
        PB.VarDepScalar("global.DEGASS", "",  "normalized DEGASS forcing"),
    
        PB.VarDepScalar("Sr_sed", "mol",  "sedimentary Sr", attributes=(:field_data=>SrIsotopeType, )),
        PB.VarDepScalar("Sr_sed_delta", "",  "d87Sr of sedimentary Sr"),
        PB.VarDepScalar("Sr_sed_norm", "",  "normalized sedimentary Sr"),    
        PB.VarContribScalar("Sr_sed_sms", "mol yr-1",  "sedimentary Sr source minus sink flux", attributes=(:field_data=>SrIsotopeType, ))
    ]

    PB.add_method_do!(
        rj, 
        do_Sr_sed,
        (PB.VarList_namedtuple(vars), ),
        p=SrIsotopeType,
    )

    return nothing
end

function PB.check_configuration(rj::ReactionSrSed, model::PB.Model)
    configok = true

    _, SrIsotopeType = PB.split_nameisotope("::SrIsotope", rj.external_parameters)
   
    if !(SrIsotopeType in (PB.ScalarData, PB.IsotopeLinear))
        @warn "ReactionSrSed.check_configuration unsupported IsotopeType $SrIsotopeType"
        configok = false
    end

    return configok
end

function do_Sr_sed(
    m::PB.ReactionMethod,
    (vars, ),
    cellrange::PB.AbstractCellRange,
    deltat
)
    rj = m.reaction
    SrIsotopeType = m.p

    # metamorphic loss from sedimentary reservoir
    if rj.pars.f_Sr_metam.v == "original"
        # mol yr-1
        Sr_metam = rj.pars.k_Sr_metam.v * vars.DEGASS[]        
    elseif rj.pars.f_Sr_metam.v == "alternative"
        # mol yr-1
        Sr_metam = rj.pars.k_Sr_metam.v * vars.DEGASS[] * vars.Sr_sed_norm[]
    else
        error("unknown f_Sr_metam $(rj.pars.f_Sr_metam.v)")
    end
    
    vars.Sr_sed_sms[] -= @PB.isotope_totaldelta(SrIsotopeType, Sr_metam, vars.Sr_sed_delta[]) 

  
    # Simplified code for time-evolution of Sr_sed due to Rb decay

    # Sediment Rb content (include a very small secular decrease in Rb abundance for consistency, ~+1.4% at 1 Ga relative to present)
    Rb_sed = PB.get_total(vars.Sr_sed[])*rj.pars.sediment_RbSr.v*exp(rj.pars.lambda_Rb.v*(-vars.tforce[]))
    # Rb decay to 87Sr
    d87Sr_sed_dt_Rb = Rb_sed*rj.pars.lambda_Rb.v
    # NB: we are (ab)using IsotopeLinear to store a linearisation of 87Sr/86Sr,
    # so use a difference to add just 87Sr ie a flux with total=0, moldelta=d87Sr_sed_dt_Rb
    Sr_sed_dt_Rb = @PB.isotope_totaldelta(SrIsotopeType, d87Sr_sed_dt_Rb, 1.0) - @PB.isotope_totaldelta(SrIsotopeType, d87Sr_sed_dt_Rb, 0.0)
    vars.Sr_sed_sms[] += Sr_sed_dt_Rb

    return nothing
end


"""
    ReactionSrLand

Calculate Sr weathering flux from land surface, given relative (normalized) basalt, granite, carbonate weathering rates and Sr isotopic composition.
"""
Base.@kwdef mutable struct ReactionSrLand{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("k_Sr_total_igw", 13e9, units="mol yr-1",
            description="total Sr weathering rate constant from igneous rocks (granite + basalt). NB: only used during configuration!!"),
        PB.ParDouble("k_Sr_basw", NaN, units="mol yr-1",
            description="Sr basalt weathering rate constant"),
        PB.ParDouble("k_Sr_granw", NaN, units="mol yr-1",
            description="Sr granite weathering rate constant from granite"),
        PB.ParString("f_Sr_sedw", "alternative", allowed_values=["original", "alternative"],
            description="functional form for sediment Sr weathering"),
        PB.ParDouble("k_Sr_sedw", 17e9, units="mol yr-1",
            description="Sr weathering rate constant from sediments"),
    )
end


function PB.register_methods!(rj::ReactionSrLand)
    _, SrIsotopeType = PB.split_nameisotope("::SrIsotope", rj.external_parameters)

    vars = [
        PB.VarDepScalar("basw_relative", "",  "Basalt weathering normalized to present"),
        PB.VarDepScalar("granw_relative", "", "Granite weathering normalized to present"),
        PB.VarDepScalar("carbw_relative", "", "Carbonate weathering normalized to present"),
    
        PB.VarContribScalar("fluxRtoOcean_Sr", "mol yr-1",  "Sr riverine flux",
            link_namestr="fluxRtoOcean.flux_Sr", attributes=(:field_data=>SrIsotopeType,)),
       
        PB.VarContribScalar("fluxLandtoSedCrust_Sr", "mol yr-1",  "Sr flux from land to sedimentary reservoirs",
            link_namestr="fluxLandtoSedCrust.flux_Sr", attributes=(:field_data=>SrIsotopeType,)),   
    
        PB.VarDepScalar("sedcrust.Sr_old_ig_delta",  "", "d87Sr old igneous rocks (granite)"),
        PB.VarDepScalar("sedcrust.Sr_new_ig_delta",  "", "d87Sr new igneous rocks (basalt)"),
        PB.VarDepScalar("sedcrust.Sr_sed_delta", "",  "d87Sr of sedimentary Sr"),
        PB.VarDepScalar("sedcrust.Sr_sed_norm", "",  "normalized sedimentary Sr"),
    ]

    PB.add_method_do!(
        rj, 
        do_Sr_land,
        (PB.VarList_namedtuple(vars), ),
        p=SrIsotopeType,
    )

    return nothing
end

function do_Sr_land(
    m::PB.ReactionMethod,
    (vars, ),
    cellrange::PB.AbstractCellRange,
    deltat
)
    rj = m.reaction
    SrIsotopeType = m.p

    # Sr weathering flux from new igneous rocks
    Sr_new_igw = @PB.isotope_totaldelta(SrIsotopeType, rj.pars.k_Sr_basw.v * vars.basw_relative[], vars.Sr_new_ig_delta[])

    # Sr weathering flux from old igneous rocks
    Sr_old_igw = @PB.isotope_totaldelta(SrIsotopeType, rj.pars.k_Sr_granw.v * vars.granw_relative[], vars.Sr_old_ig_delta[])

    # Sr weathering from sedimentary reservoir
    if rj.pars.f_Sr_sedw.v == "original"
        # mol yr-1
        Sr_sedw_tot = rj.pars.k_Sr_sedw.v * vars.carbw_relative[]        
    elseif rj.pars.f_Sr_sedw.v == "alternative"
        # mol yr-1
        Sr_sedw_tot = rj.pars.k_Sr_sedw.v * vars.carbw_relative[] * vars.Sr_sed_norm[]
    else
        error("unknown f_Sr_sedw $(rj.pars.f_Sr_sedw.v)")
    end    
    Sr_sedw = @PB.isotope_totaldelta(SrIsotopeType, Sr_sedw_tot, vars.Sr_sed_delta[]) 

    # Riverine flux
    vars.fluxRtoOcean_Sr[]        += Sr_new_igw + Sr_old_igw + Sr_sedw 

    # Flux to sedimentary Sr reservoir (-ve)
    vars.fluxLandtoSedCrust_Sr[]  -= Sr_sedw

    return nothing
end


"""
    ReactionSrOceanfloor

Calculate evolution of an ocean`Sr` reservoir due to ocean burial (assumed proportional to carbonate burial), 
seafloor weathering, and hydrothermal (mantle) input fluxes.
"""
Base.@kwdef mutable struct ReactionSrOceanfloor{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("k_Sr_mantle", 7.3e9, units="mol yr-1",
            description="Sr mantle (hydrothermal) input"),
        PB.ParDouble("k_Sr_sfw", NaN, units="mol yr-1",
            description="Sr seafloor weathering output"),
    
        PB.ParString("f_Sr_sedb", "carbburial", allowed_values=["carbburial", "silwcarbw"],
            description="functional form for Sr burial"),
        PB.ParDouble("k_Sr_sedb", NaN, units="mol yr-1",
            description="Sr burial output"),
        PB.ParDouble("k_mccb_0", NaN, units="mol yr-1",
            description="Carbonate burial rate to normalize Sr burial output for f_Sr_sedb=\"carbburial\""),
    )
end


function PB.register_methods!(rj::ReactionSrOceanfloor)
    _, SrIsotopeType = PB.split_nameisotope("::SrIsotope", rj.external_parameters)

    vars = [
        PB.VarDepScalar("global.DEGASS",        "",         "normalized DEGASS forcing"),
        PB.VarDepScalar("land.silwcarbw_relative","",       "normalized land silicate+carbonate weathering"),
        PB.VarDepScalar("sfw_relative",         "",         "normalized seafloor weathering"),
        
        PB.VarContrib("solutefluxOceanfloor_Sr","mol yr-1", "Sr oceanfloor solute flux",
            link_namestr="fluxOceanfloor.soluteflux_Sr", attributes=(:field_data=>SrIsotopeType,)),
       
        PB.VarContrib("fluxOceanBurial_Sr",     "mol yr-1", "Sr ocean burial flux",
            link_namestr="fluxOceanBurial.flux_Sr", attributes=(:field_data=>SrIsotopeType,)),
        PB.VarDep("fluxOceanBurial_Ccarb",      "mol yr-1", "carbonate ocean burial flux",
            link_namestr="fluxOceanBurial.flux_Ccarb"),
           
        PB.VarDepScalar("ocean.Sr_norm",        "",         "normalized ocean Sr"),
        PB.VarDep("ocean.oceanfloor.Sr_delta",  "",         "ocean d87Sr"), 
        PB.VarDep("sedcrust.Sr_mantle_delta",   "",         "mantle d87Sr"), 
    ]

    PB.add_method_do!(
        rj, 
        do_Sr_oceanfloor,
        (PB.VarList_namedtuple(vars), ),
        p=SrIsotopeType,
    )

    return nothing
end

function do_Sr_oceanfloor(
    m::PB.ReactionMethod,
    (vars, ),
    cellrange::PB.AbstractCellRange,
    deltat
)
    rj = m.reaction
    SrIsotopeType = m.p

    r_nfloorcells = 1.0/PB.get_length(rj.domain) # fraction of flux for each oceanfloor cell

    # Sr mantle (hydrothermal) input  
    Sr_mantle_total = @PB.isotope_totaldelta(SrIsotopeType, rj.pars.k_Sr_mantle.v * vars.DEGASS[], vars.Sr_mantle_delta[])
    @inbounds for i in cellrange.indices
        vars.solutefluxOceanfloor_Sr[i] += r_nfloorcells*Sr_mantle_total
    end

    # Sr seafloor weathering output
    Sr_sfw_total = rj.pars.k_Sr_sfw.v * vars.sfw_relative[] * vars.Sr_norm[]
    # @Infiltrator.infiltrate
    @inbounds for i in cellrange.indices
        vars.solutefluxOceanfloor_Sr[i] -= r_nfloorcells*@PB.isotope_totaldelta(SrIsotopeType, Sr_sfw_total, vars.Sr_delta[i])
    end

    # Sr burial
    if rj.pars.f_Sr_sedb.v == "silwcarbw"
        # COPSE formulation: assumes mccb ~ silw + carbw
        # mol yr-1
        Sr_sedb_total = rj.pars.k_Sr_sedb.v * vars.silwcarbw_relative[] * vars.Sr_norm[]
        @inbounds for i in cellrange.indices
            Sr_sedb = r_nfloorcells*@PB.isotope_totaldelta(SrIsotopeType, Sr_sedb_total, vars.Sr_delta[i])
            vars.solutefluxOceanfloor_Sr[i] -= Sr_sedb
            vars.fluxOceanBurial_Sr[i]      += Sr_sedb
        end   
    elseif rj.pars.f_Sr_sedb.v == "carbburial"
        # Proportional to actual carbonate burial rate
        # mol mol-1
        SrCarbRatio = rj.pars.k_Sr_sedb.v/rj.pars.k_mccb_0.v * vars.Sr_norm[]
        @inbounds for i in cellrange.indices
            Sr_sedb = r_nfloorcells*@PB.isotope_totaldelta(SrIsotopeType, rj.fluxOceanBurial_Ccarb[i]*SrCarbRatio, vars.Sr_delta[i])
            vars.solutefluxOceanfloor_Sr[i] -= Sr_sedb
            vars.fluxOceanBurial_Sr[i]      += Sr_sedb
        end   
       
    else
        error("unknown f_Sr_sedb $(rj.pars.f_Sr_sedb.v)")
    end

    return nothing
end


"""
    set_Sr_fluxes_steady_state!(
        rct_Sr_land::ReactionSrLand,
        rct_Sr_oceanfloor::ReactionSrOceanfloor,
        basw, granw, carbw, sfw
    )

Set Sr weathering fluxes for steady-state, given present-day basalt, granite, carbonate and seafloor weathering rates (mol yr-1)
"""
function set_Sr_fluxes_steady_state!(
    rct_Sr_land::ReactionSrLand,
    rct_Sr_oceanfloor::ReactionSrOceanfloor,
    basw, granw, carbw, sfw
)

    # Partition Sr weathering from igneous rocks 
    silw = basw + granw
    PB.setvalue!(rct_Sr_land.pars.k_Sr_basw,  basw/silw*rct_Sr_land.pars.k_Sr_total_igw.v)
    PB.setvalue!(rct_Sr_land.pars.k_Sr_granw,  granw/silw*rct_Sr_land.pars.k_Sr_total_igw.v)

    # Total Sr inputs
    Sr_total_inputs = rct_Sr_land.pars.k_Sr_sedw.v + rct_Sr_land.pars.k_Sr_total_igw.v + rct_Sr_oceanfloor.pars.k_Sr_mantle.v
    @info "set_Sr_fluxes_steady_state!:  Sr_total_inputs = $Sr_total_inputs (mol Sr yr-1)"

    # Set total present-day carbonate burial (for normalisation of Sr burial flux) TODO omitting seafloor weathering?
    PB.setvalue!(rct_Sr_oceanfloor.pars.k_mccb_0,  silw + carbw)

    # Set Sr seafloor weathering output assuming fraction of Sr to sfw is same as fraction of carbonate
    PB.setvalue!(rct_Sr_oceanfloor.pars.k_Sr_sfw, Sr_total_inputs * sfw / ( silw + carbw + sfw))
    @info "set_Sr_fluxes_steady_state!:  seafloor weathering k_Sr_sfw = $(rct_Sr_oceanfloor.pars.k_Sr_sfw.v) (mol Sr yr-1)"

    # Set Sr burial output for steady state
    PB.setvalue!(rct_Sr_oceanfloor.pars.k_Sr_sedb, Sr_total_inputs - rct_Sr_oceanfloor.pars.k_Sr_sfw.v)
    @info "set_Sr_fluxes_steady_state!:  ocean burial k_Sr_sedb = $(rct_Sr_oceanfloor.pars.k_Sr_sedb.v) (mol Sr yr-1)"


    return nothing
end


"Install create_reactionXXX factories when module imported"
function __init__()    
    PB.add_reaction_factory(ReactionSrMantleCrust)
    PB.add_reaction_factory(ReactionSrSed)
    PB.add_reaction_factory(ReactionSrLand)
    PB.add_reaction_factory(ReactionSrOceanfloor)

    return nothing
end

end
