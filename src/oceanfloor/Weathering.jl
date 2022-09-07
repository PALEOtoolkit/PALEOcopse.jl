# -*- coding: utf-8 -*-
module Weathering


import PALEOboxes as PB
using PALEOboxes.DocStrings

"""
    ReactionSeafloorWeathering

COPSE Reloaded (Lenton etal 2018) seafloor weathering (basalt carbonation).

In a spatially resolved oceanfloor Domain, the total flux is calculated in the same way as for the 0D COPSE model,
and is then distributed over a range of oceanfloor cells according the Paramter `sfw_distribution_method`.

Fluxes are added to flux couplers:
- `fluxOceanfloor`: ocean solute fluxes
- `fluxOceanBurial`: carbonate burial flux

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionSeafloorWeathering{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParString("f_sfw_opt",     "sfw_temp",  allowed_values=["sfw_temp", "sfw_strong", "sfw_noT", "mills2014pCO2"],
            description="functional form for seafloor weathering rate temperature dependence"),
        PB.ParDouble("f_sfw_alpha",   0.23,     units="",
            description="power law for f_sfw_opt=mills2014pCO2 pCO2 dependence: wildly uncertain"),
        PB.ParString("f_sfw_force",   "DEGASS",  allowed_values=["None", "DEGASS"],
            description="seafloor weathering rate tectonic forcing"),
        
        PB.ParDouble("k_sfw",   3e12,       units="mol/yr",
            description="seafloor weathering rate"),
        PB.ParString("sfw_distribution_method", "Parameter", allowed_values=["Parameter", "Depth"],
            description="method to set sfw distribution: "*
                        "Parameter - use sfw_distribution parameters, "*
                        "Depth - distribute proportional to area within specified depth range"),
        PB.ParDoubleVec("sfw_distribution", [1.0],
            description="per box distribution of seafloor weathering (length=Oceanfloor Domain size, must sum to 1.0)"),
        PB.ParDoubleVec("sfw_depthrange", [-3000.0, -1000.0], units="m",
            description="depth range for distribution of seafloor weathering (lower, upper)" ),

        PB.ParString("f_sfw_d13C",   "delta_DIC",  allowed_values=["delta_DIC", "delta_mccb"],
            description="d13C delta"),

        # Isotopes
        PB.ParType(PB.AbstractData, "CIsotope", PB.ScalarData,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable carbon isotopes and specify isotope type"),
    )
end


function PB.register_methods!(rj::ReactionSeafloorWeathering)

    CIsotopeType = rj.pars.CIsotope[]
    PB.setfrozen!(rj.pars.CIsotope)
    vars = [
        PB.VarDepScalar("global.RHOSFW",          "", "seafloor weathering-specific additional forcing (usually 1.0)"),
        PB.VarDepScalar("(global.TEMP)",          "K", "global mean temperature"),
        PB.VarDepScalar("(global.DEGASS)",      "", "degass forcing"),
        PB.VarDepScalar("(atm.pCO2PAL)",          "PAL", "atmospheric pCO2PAL"),
        PB.VarDep("(ocean.oceanfloor.DIC_delta)",          "per mil", "d13C DIC"),
        PB.VarDepScalar("(ocean.D_mccb_DIC)",   "per mil",  "d13C marine calcite burial relative to ocean DIC"),

        PB.VarPropScalar("sfw_relative",  "", "normalized seafloor weathering"),
        PB.VarProp("sfw",     "mol yr-1", "seafloor weathering flux",
            attributes=(:field_data=>CIsotopeType, :calc_total=>true,)),
        # PB.VarPropScalar("sfw_total",     "mol yr-1", "total seafloor weathering flux"
    ]        

    
    # flux couplers
    fluxOceanBurial = PB.Fluxes.FluxContrib(
        "fluxOceanBurial.flux_",
        ["Ccarb::$CIsotopeType"],
        isotope_data=Dict())
    
    fluxOceanfloorSolute = PB.Fluxes.FluxContrib(
        "fluxOceanfloor.soluteflux_",
        ["DIC::$CIsotopeType"],
        isotope_data=Dict())
    
    if rj.pars.sfw_distribution_method[] == "Depth"
        grid_vars = [
            PB.VarDep("ocean.oceanfloor.zlower",    "m",    "oceanfloor depth"),
            PB.VarDep("oceanfloor.Afloor",          "m^2",  "horizontal area of seafloor at base of box"),
        ]
        PB.add_method_setup!(
            rj, 
            set_distributionfromdepths,
            (PB.VarList_namedtuple(grid_vars),)
        )
    end

    PB.add_method_do!(
        rj, 
        do_seafloor_weathering,
        (
            PB.VarList_namedtuple_fields(fluxOceanBurial), 
            PB.VarList_namedtuple_fields(fluxOceanfloorSolute),
            PB.VarList_namedtuple(vars),
        ),
        p=CIsotopeType,
    )
 
    PB.add_method_do_totals_default!(rj)

    PB.add_method_initialize_zero_vars_default!(rj) # for total variable(s)

    return nothing
end

function PB.check_configuration(rj::ReactionSeafloorWeathering, model::PB.Model)
    configok = true
    
    if rj.pars.sfw_distribution_method[] == "Parameter"
        PB.check_parameter_sum(rj.pars.sfw_distribution, PB.get_length(rj.domain)) || (configok = false)
    elseif rj.pars.sfw_distribution_method[] == "Depth"
        (length(rj.pars.sfw_depthrange) == 2 && rj.pars.sfw_depthrange[1] < rj.pars.sfw_depthrange[2]) || 
            (@warn "invalid sfw_depthrange $(rj.pars.sfw_depthrange.v)"; configok = false)
    else
        @warn("unknown sfw_distribution_method $(rj.pars.sfw_distribution_method[])")
    end

    return configok
end


function do_seafloor_weathering(
    m::PB.ReactionMethod,
    pars,
    (fluxOceanBurial, fluxOceanfloorSolute, vars), 
    cellrange::PB.AbstractCellRange,
    deltat
)
    CIsotopeType = m.p

    # Tectonic dependence
    if pars.f_sfw_force[]  == "None"
        force_sfw = vars.RHOSFW[]
    elseif pars.f_sfw_force[] == "DEGASS"
        force_sfw = vars.RHOSFW[] * vars.DEGASS[]
    else
        error("unknown f_sfw_force ", pars.f_sfwforce[])
    end

    # Temperature / CO2 dependence
    # NB: "sfw_Tbotw" and "sfw_temp" temperature-based functions are _very_ similar 
    if      pars.f_sfw_opt[] == "mills2014pCO2"       # Mills (2014) PNAS 10.1073/pnas.1321679111
        f_sfw = vars.pCO2PAL[]^pars.f_sfw_alpha[]        
    elseif  pars.f_sfw_opt[] == "sfw_Tbotw"           # Stuart sfw function for ox weath model
        # normalised to bottom-water temperature (global TEMP(in C) - 12.5)
        # from S&G high-lat T = 2.5C for Temp = 15C, don't go below freezing        
        TbotwC = max(vars.TEMP[] - PB.Constants.k_CtoK - 12.5, 0) 
        f_sfw = exp(0.066*(TbotwC - 2.5)) # activation energy 41 kJ mol^{-1}
        # NB: oxweath had normalisation error, gave f_sfw = exp(0.165) = 1.18 for present day TEMP = 15C        
    elseif  pars.f_sfw_opt[] == "sfw_temp"            # Josh Improved sfw function considering temperature
        f_sfw =  exp(0.0608*(vars.TEMP[] - 288)) # 42KJ/mol activation energy assumed as with terrestrial basalt      
    elseif  pars.f_sfw_opt[] == "sfw_strong"          # Coogan&Dosso 92 kJ/mol apparent activation energy
        f_sfw =  exp(0.1332*(vars.TEMP[] - 288.15))        
    elseif  pars.f_sfw_opt[] == "sfw_noT"             # No temperature dependence
        f_sfw =  1.0        
    else
        error("unrecognized f_sfw_opt ", pars.f_sfw_opt[])
    end
            
    # calculate relative rate of seafloor weathering 
    vars.sfw_relative[] =   force_sfw * f_sfw

    # total flux, without isotope contribution
    sfw_total_noisotope = vars.sfw_relative[] * pars.k_sfw[] 

    # isotope composition from local DIC with optional offset 
    if CIsotopeType <: PB.AbstractIsotopeScalar
        if      pars.f_sfw_d13C[] == "delta_DIC"
            d13C_offset = 0.0
        elseif  pars.f_sfw_d13C[] == "delta_mccb"
            d13C_offset = vars.D_mccb_DIC[]
        else
            error("unknown f_sfw_d13C ", pars.f_sfw_d13C[])
        end
    end

    # distribute total flux to seafloor boxes, adding isotope composition
    @inbounds for i in cellrange.indices
        vars.sfw[i] = @PB.isotope_totaldelta(
            CIsotopeType, 
            sfw_total_noisotope*pars.sfw_distribution[i], 
            vars.DIC_delta[i] + d13C_offset
        )
        fluxOceanBurial.Ccarb[i]    += vars.sfw[i]
        fluxOceanfloorSolute.DIC[i] -= vars.sfw[i]
    end

    return nothing
end


"set distribution of sfw flux for equal flux per unit area for boxes in specified depth range"
function set_distributionfromdepths(
    m::PB.ReactionMethod,
    pars,
    (vars, ), 
    cellrange::PB.AbstractCellRange, 
    attribute_name
)
    rj = m.reaction
    attribute_name == :setup || return nothing

    @info "$(PB.fullname(m)) ReactionSeafloorWeathering:"

    PB.get_length(rj.domain) == length(cellrange.indices) || 
        error("ReactionSeafloorWeathering $(PB.fullname(rj)) set_distributionfromdepths cellrange does not cover whole Domain")

    sfw_distribution = zeros(PB.get_length(rj.domain))
    depth_lower, depth_upper = pars.sfw_depthrange

    nsfw = 0.0
    for i in eachindex(sfw_distribution)      
        zlower = vars.zlower[i]
        if zlower >= depth_lower && zlower <= depth_upper
            sfw_distribution[i] = vars.Afloor[i]
            nsfw += 1
        end
    end
    
    normsum = sum(sfw_distribution)
    Atot = sum(vars.Afloor)
    @info "    sfw distributed among $nsfw of $(length(sfw_distribution)) boxes, area $normsum m^2 of $Atot m^2"
    PB.setvalue!(pars.sfw_distribution, sfw_distribution./normsum)

    return nothing
end


end # module
