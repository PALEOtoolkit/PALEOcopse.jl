module COPSEForcings


import MAT   # Matlab file access
import Interpolations

import XLSX  # Excel file access
import DataFrames 

import PALEOboxes as PB
using PALEOboxes.DocStrings

import PALEOcopse


"""
    ReactionForce_CK_Solar

`SOLAR` forcing (incident solar flux at Earth, W m^-2)

# Parameters
$(PARS)

# Methods and Variables
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionForce_CK_Solar <:  PB.AbstractReaction
    base::PB.ReactionBase

    solarpresentday::Float64 = 1368.0     # W m^-2 present-day insolation
end

function PB.register_methods!(rj::ReactionForce_CK_Solar)
    
    vars = [
        PB.VarDepScalar("tforce", "yr",  "historical time at which to apply forcings, present = 0 yr"),
        PB.VarPropScalar("SOLAR", "W m-2",  "incident solar flux at Earth")
    ]

    PB.add_method_do!(rj, do_force_CK_solar, (PB.VarList_namedtuple(vars), ))
    return nothing
end

function do_force_CK_solar(m::PB.ReactionMethod, (vars, ), cellrange::PB.AbstractCellRange, deltat)
    rj = m.reaction
    vars.SOLAR[] =  rj.solarpresentday / (1.0 + 0.38*(-vars.tforce[]/4.55e9))
    return nothing
end


"""
    ReactionForce_UDWEbergman2004

COPSE Bergman(2004) forcings from file.
    
Provides:
- `UPLIFT`  tectonic uplift. GEOCARB II (Berner 1994) based on Sr isotopes.
- `DEGASS`  metamorphic and volcanic degassing.
       Bergman etal (2004), same as GEOCARB II (Berner 1994):
       Engerbretson etal (1992) seafloor subduction (spreading) rate 0 - 150Ma,
       Gaffin (1987) Paleo-sealevel based 570 - 150Ma
- `EVO`  land plant evolution and colonisation.
      Bergman etal (2004)
- `W`  biological enhancement of weathering.
     Bergman etal (2004)

Interpolates and optionally applies naive extrapolation of forcings into (constant) Precambrian and future

# Methods and Variables
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionForce_UDWEbergman2004 <:  PB.AbstractReaction
    base::PB.ReactionBase

    C_514_forcings = Dict{String, Vector{Float64}}()

    interp_UPLIFT = nothing
    interp_DEGASS = nothing
    interp_EVO    = nothing
    interp_W      = nothing
end


function PB.register_methods!(rj::ReactionForce_UDWEbergman2004)
    forcingfile = joinpath(PALEOcopse.Forcings.srcdir(), "copse_forcings.mat")
    @info "ReactionForce_UDWEbergman2004: loading U,D,W,E forcings from COPSE datafile $forcingfile"
    
    vars = MAT.matread(forcingfile)
    copse_forcings = vars["copse_forcings"]

    # convert time to forwards in years
    rj.C_514_forcings["Tyr"] = -copse_forcings[:, 1]
    rj.C_514_forcings["UPLIFT"] = copse_forcings[:, 2]
    rj.C_514_forcings["DEGASS"] = copse_forcings[:, 3]
    rj.C_514_forcings["W"] = copse_forcings[:, 4]
    rj.C_514_forcings["EVO"] = copse_forcings[:, 5]

    # create interpolation objects
    rj.interp_UPLIFT = Interpolations.linear_interpolation(
        rj.C_514_forcings["Tyr"], rj.C_514_forcings["UPLIFT"], extrapolation_bc = 1.0
    )
    rj.interp_DEGASS = Interpolations.linear_interpolation(
        rj.C_514_forcings["Tyr"], rj.C_514_forcings["DEGASS"], extrapolation_bc = 1.0
    )
    rj.interp_EVO = Interpolations.linear_interpolation(
        rj.C_514_forcings["Tyr"], rj.C_514_forcings["EVO"], extrapolation_bc = Interpolations.Flat()
    )
    rj.interp_W = Interpolations.linear_interpolation(
        rj.C_514_forcings["Tyr"], rj.C_514_forcings["W"], extrapolation_bc = Interpolations.Flat()
    )

    vars = [
        PB.VarDepScalar("tforce",     "yr",  "historical time at which to apply forcings, present = 0 yr"),
        PB.VarPropScalar("UPLIFT",    "",  "tectonic uplift normalized to present. GEOCARB II (Berner 1994) based on Sr isotopes."),
        PB.VarPropScalar("DEGASS",    "",  "metamorphic and volcanic degassing normalized to present"),
        PB.VarPropScalar("EVO",       "",  "land plant evolution and colonisation"),
        PB.VarPropScalar("W",         "",  "biological enhancement of weathering"),
    ]

    PB.add_method_do!(
        rj,
        do_force_UDWEbergman2004,
        (PB.VarList_namedtuple(vars), ),
        p = (UPLIFT=rj.interp_UPLIFT, DEGASS=rj.interp_DEGASS, EVO=rj.interp_EVO, W=rj.interp_W),
    )

    return nothing
end

function do_force_UDWEbergman2004(m::PB.ReactionMethod, (vars, ), cellrange::PB.AbstractCellRange, deltat)

    interp = m.p

    tforce = vars.tforce[]

    vars.UPLIFT[] =  interp.UPLIFT(tforce)
    vars.DEGASS[] =  interp.DEGASS(tforce)
    vars.W[]      =  interp.W(tforce)
    vars.EVO[]    =  interp.EVO(tforce)

    return nothing
end


"""
    ReactionForce_Bbergman2004

`Bforcing` forcing, calcerous plankton evolution (COPSE 2004)

# Methods and Variables
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionForce_Bbergman2004 <:  PB.AbstractReaction
    base::PB.ReactionBase

    Btimes =    [-150e6,      -140e6,      -110e6,      -90e6,       -50e6,   	-10e6]
    Bvals  =    [0.75,       0.83776596, 0.90990691,  0.96110372,	0.98902926,	 1.0]

    interp_B = Interpolations.linear_interpolation(
        Btimes,
        Bvals,
        extrapolation_bc=Interpolations.Flat()
    )
end

function PB.register_methods!(rj::ReactionForce_Bbergman2004)

    vars = [
        PB.VarDepScalar("tforce", "yr",  "historical time at which to apply forcings, present = 0 yr"),
        PB.VarPropScalar("Bforcing", "",  "calcerous plankton evolution")
    ]

    PB.add_method_do!(rj, do_force_Bbergman2004, (PB.VarList_namedtuple(vars), ), p=rj.interp_B)

    return nothing
end


function do_force_Bbergman2004(m::PB.ReactionMethod, (vars, ), cellrange::PB.AbstractCellRange, deltat)

    interp_B = m.p
    vars.Bforcing[] =  interp_B(vars.tforce[])
    return nothing
end


"""
    ReactionForce_CPlandrelbergman2004

`CPland_relative` forcing.
CP land burial ratio doubles in permo-carboniferous (COPSE)

# Methods and Variables
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionForce_CPlandrelbergman2004 <: PB.AbstractReaction
    base::PB.ReactionBase
 
    CPlandtimes =    [-355e6,      -345e6,      -290e6,      -280e6]
    CPlandvals  =    [1.0,          2.0,        2.0,        1.0]

    interp_CPland = Interpolations.linear_interpolation(
        CPlandtimes, 
        CPlandvals, 
        extrapolation_bc = Interpolations.Flat()
    )
end

function PB.register_methods!(rj::ReactionForce_CPlandrelbergman2004)

    vars = [
        PB.VarDepScalar("tforce", "yr",  "historical time at which to apply forcings, present = 0 yr"),
        PB.VarPropScalar("CPland_relative", "",  "normalized land CP burial ratio"),
    ]

    PB.add_method_do!(
        rj,
        do_force_CPlandrelbergman2004,
        (PB.VarList_namedtuple(vars), ),
        p=rj.interp_CPland
    )
end

function do_force_CPlandrelbergman2004(m::PB.ReactionMethod, (vars, ), cellrange::PB.AbstractCellRange, deltat)
    interp_CPland = m.p

    vars.CPland_relative[] =  interp_CPland(vars.tforce[])

    return nothing
end

const LINEARINTERPOLATION_TEMPLATE = Interpolations.linear_interpolation(
    [0.0, 1.0], 
    [NaN, NaN],
    extrapolation_bc = Interpolations.Flat()
)

"""
    ReactionForce_spreadsheet

Generic forcing interpolated from values in spreadsheet.

Spreadsheet should contain a table of numeric data in sheet `sheetname`,
with a single header row, and time in `timecolumn` and data in `datacolumn`.

Time from `timecolumn` is multiplied by `timemultiplier` to convert to PALEO model time.

Linearly interpolates in time `tforce` to set Variable `forcename`, and extrapolates into past and future, either to values
in `extrap_value_past`, `extrap_value_future`, or to the closest spreadsheet value if these Parameters are `NaN`.

# Implementation
Spreadsheet data is read using the Julia XLSX and DataFrames packages with

    xf = XLSX.readxlsx(forcingfile)
    df = XLSX.eachtablerow(xf[sheetname]) |> DataFrames.DataFrame

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionForce_spreadsheet{P} <:  PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParString("forcename", "FORCE",
            description="name of forcing (this will be a variable name in global Domain)"),
        PB.ParString("datafolder", PALEOcopse.Forcings.srcdir(),
            description="folder for spreadsheet with forcing data"),
        PB.ParString("datafile", "", 
            description="spreadsheet with forcing data (path relative to 'datafolder')"),
        PB.ParString("sheetname", "Sheet1",
            description="sheet name in spreadsheet"),
        PB.ParInt("timecolumn", 1,
            description="column with time data"),
        PB.ParDouble("timemultiplier", -1.0e6,
            description="factor to multiply 'timecolumn' by to convert to yr after present day (so times in past are -ve)"),
        PB.ParInt("datacolumn", 2,
            description="column with forcing data"),
        PB.ParDouble("extrap_value_past", 1.0, units="", 
            description="extrapolate value for tforce before earliest value in spreadsheet (NaN to use earliest value)"),
        PB.ParDouble("extrap_value_future", 1.0, units="", 
            description="extrapolate value for tforce after latest value in spreadsheet (NaN to use latest value"),
    )

    forcing_data    = DataFrames.DataFrame()  # raw data from spreadsheet

    force_times::Vector{Float64} = Float64[] # times used for interpolation
    force_values::Vector{Float64} = Float64[] # values used for interpolation
    interp_FORCE::typeof(LINEARINTERPOLATION_TEMPLATE) = LINEARINTERPOLATION_TEMPLATE
    
end

function PB.register_methods!(rj::ReactionForce_spreadsheet)
  
    var_tforce = PB.VarDepScalar("tforce",     "yr",  
        "historical time at which to apply forcings, present = 0 yr")
    var_FORCE = PB.VarPropScalar(rj.pars.forcename[],    "",  
        "forcing interpolated from spreadsheet")
    PB.setfrozen!(rj.pars.forcename)

    PB.add_method_setup!(
        rj, 
        setup_force_spreadsheet,
        (),
    )

    PB.add_method_do!(
        rj, 
        do_force_spreadsheet,
        (PB.VarList_single(var_tforce), PB.VarList_single(var_FORCE)),
    )

    return nothing
end

function setup_force_spreadsheet(m::PB.ReactionMethod, pars, (), cellrange::PB.AbstractCellRange, attribute_name)
    rj = m.reaction

    attribute_name == :setup || return nothing

    forcingfile = joinpath(pars.datafolder[], pars.datafile[])

    io = IOBuffer()
    println(io, "setup_force_spreadsheet ReactionForce_spreadsheet $(PB.fullname(rj)): ")
    println(io, "    loading $(pars.forcename[]) forcing from 'datafolder/datafile'='$(forcingfile)'")
    rj.forcing_data = _read_xlsx(io, forcingfile, sheetname=pars.sheetname[])

    sp_times = pars.timemultiplier[]*rj.forcing_data[:, pars.timecolumn[]]
    println(io, "    'tforce' from $(pars.timemultiplier[]) * column $(pars.timecolumn[]) ($(names(rj.forcing_data)[pars.timecolumn[]]))")

    sp_values = Float64.(rj.forcing_data[:, pars.datacolumn[]])
    println(io, "    '$(pars.forcename[])' from column $(pars.datacolumn[]) ($(names(rj.forcing_data)[pars.datacolumn[]]))")

    # sort in ascending time order
    sp_perm = sortperm(sp_times)
    rj.force_times = sp_times[sp_perm]
    rj.force_values = sp_values[sp_perm]

    extrap_past = isnan(pars.extrap_value_past[]) ? "earlist value in spreadsheet = $(first(rj.force_values))" : "'extrap_value_past' = $(pars.extrap_value_past[])"
    println(io, "    extrapolating out-of-range tforce < $(first(rj.force_times)) (yr) to $extrap_past")
    extrap_future = isnan(pars.extrap_value_future[]) ? "latest value in spreadsheet = $(last(rj.force_values))" : "'extrap_value_future' = $(pars.extrap_value_future[])"
    println(io, "    extrapolating out-of-range tforce > $(last(rj.force_times)) (yr) to $extrap_future")
   
    # create interpolation object
    rj.interp_FORCE = Interpolations.linear_interpolation(
        rj.force_times, rj.force_values, 
        extrapolation_bc = Interpolations.Flat() # only used for extrap_value_past, future == NaN
    )

    @info String(take!(io))

    return nothing
end

function do_force_spreadsheet(m::PB.ReactionMethod, pars, (var_tforce, var_FORCE), cellrange::PB.AbstractCellRange, deltat)
    rj = m.reaction

    tforce = var_tforce[]

    if tforce < first(rj.force_times) && !isnan(pars.extrap_value_past[])
        var_FORCE[] = pars.extrap_value_past[]
    elseif tforce > last(rj.force_times) && !isnan(pars.extrap_value_future[])
        var_FORCE[] = pars.extrap_value_future[]
    else
        # extrapolation_bc = Flat() will extrapolate to first/last constant value
        var_FORCE[] =  rj.interp_FORCE(tforce)
    end
    
    return nothing
end

function _read_xlsx(io, forcingfile; sheetname="Sheet1")

    println(io, "    read_xlsx: spreadsheet $(forcingfile) sheet $(sheetname)")

    xf = XLSX.readxlsx(forcingfile)

    # Read using built-in iterator - assumes single header row, automatic type conversion
    df = XLSX.eachtablerow(xf[sheetname]) |> DataFrames.DataFrame

    println(io, "    read_xlsx read $(DataFrames.describe(df))")

    return df
end

end # module
