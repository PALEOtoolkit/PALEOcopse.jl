module COPSEForcings


import MAT   # Matlab file access
import Interpolations

import XLSX  # Excel file access
import DataFrames 

import PALEOboxes as PB
import PALEOcopse


"""
    ReactionForce_CK_Solar

`SOLAR` forcing (incident solar flux at Earth, W m^-2)
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
    rj.interp_UPLIFT = Interpolations.LinearInterpolation(
        rj.C_514_forcings["Tyr"], rj.C_514_forcings["UPLIFT"], extrapolation_bc = 1.0
    )
    rj.interp_DEGASS = Interpolations.LinearInterpolation(
        rj.C_514_forcings["Tyr"], rj.C_514_forcings["DEGASS"], extrapolation_bc = 1.0
    )
    rj.interp_EVO = Interpolations.LinearInterpolation(
        rj.C_514_forcings["Tyr"], rj.C_514_forcings["EVO"], extrapolation_bc = Interpolations.Flat()
    )
    rj.interp_W = Interpolations.LinearInterpolation(
        rj.C_514_forcings["Tyr"], rj.C_514_forcings["W"], extrapolation_bc = Interpolations.Flat()
    )

    vars = [
        PB.VarDepScalar("tforce",     "yr",  "historical time at which to apply forcings, present = 0 yr"),
        PB.VarPropScalar("UPLIFT",    "",  "tectonic uplift normalized to present"),
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
"""
Base.@kwdef mutable struct ReactionForce_Bbergman2004 <:  PB.AbstractReaction
    base::PB.ReactionBase

    Btimes =    [-150e6,      -140e6,      -110e6,      -90e6,       -50e6,   	-10e6]
    Bvals  =    [0.75,       0.83776596, 0.90990691,  0.96110372,	0.98902926,	 1.0]

    interp_B = Interpolations.LinearInterpolation(
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
"""
Base.@kwdef mutable struct ReactionForce_CPlandrelbergman2004 <: PB.AbstractReaction
    base::PB.ReactionBase
 
    CPlandtimes =    [-355e6,      -345e6,      -290e6,      -280e6]
    CPlandvals  =    [1.0,          2.0,        2.0,        1.0]

    interp_CPland = Interpolations.LinearInterpolation(
        CPlandtimes, 
        CPlandvals, 
        extrapolation_bc = Interpolations.Flat()
    )
end

function PB.register_methods!(rj::ReactionForce_CPlandrelbergman2004)

    vars = [
        PB.VarDepScalar("tforce", "yr",  "historical time at which to apply forcings, present = 0 yr"),
        PB.VarProp("CPland_relative", "",  "normalized land CP burial ratio"),
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


"""
    ReactionForce_spreadsheet

Generic forcing interpolated from values in spreadsheet. Spreadsheet should contain a table of numeric data in Sheet1,
with a single header row, and time in Ma as the first column.

Interpolates and optionally applies naive extrapolation of forcings into (constant) Precambrian and future
"""
Base.@kwdef mutable struct ReactionForce_spreadsheet{P} <:  PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParString("forcename", "FORCE",
            description="name of forcing (this will be a variable name in global Domain)"),
        PB.ParString("datafile", "", 
            description="spreadsheet with forcing data (path relative to $(PALEOcopse.Forcings.srcdir()))"),
        PB.ParInt("timecolumn", 1,
            description="column with time data"),
        PB.ParInt("datacolumn", 2,
            description="column with forcing data"),
        PB.ParDouble("extrap_value", 1.0, units="", 
            description="extrapolate value for out-of-range tforce"),
    )

    forcing_data    = DataFrames.DataFrame()  # raw data from spreadsheet

    interp_FORCE = nothing
    
end

function read_xlsx(forcingfile; sheetname="Sheet1")

    @info "read_xlsx: spreadsheet $(forcingfile) sheet $(sheetname)"

    xf = XLSX.readxlsx(forcingfile)
    # read data assuming first row is column headers
    # data, column_labels = XLSX.gettable(xf[sheetname])
    # convert to float64 and create DataFrame
    # data_float64 = [Float64.(d) for d in data]
    # df = DataFrames.DataFrame(data_float64, column_labels)

    # Read using built-in iterator - assumes single header row, automatic type conversion
    df = XLSX.eachtablerow(xf[sheetname]) |> DataFrames.DataFrame

    @info "read_xlsx read $(DataFrames.describe(df))"

    return df
end

function PB.register_methods!(rj::ReactionForce_spreadsheet)
  
    forcingfile = joinpath(PALEOcopse.Forcings.srcdir(), rj.pars.datafile.v)
    @info "ReactionForce_spreadsheet: $(PB.fullname(rj)) loading $(rj.pars.forcename.v) forcing from datafile $(forcingfile)"
    rj.forcing_data = read_xlsx(forcingfile)
    @info "ReactionForce_spreadsheet:  $(PB.fullname(rj)) $(rj.pars.forcename.v) from column $(rj.pars.datacolumn.v) ($(names(rj.forcing_data)[rj.pars.datacolumn.v]))"
    @info "ReactionForce_spreadsheet:  $(PB.fullname(rj)) extrapolating out-of-range tforce to $(rj.pars.forcename.v) = $(rj.pars.extrap_value.v)"
    
    # create interpolation object
    rj.interp_FORCE = Interpolations.LinearInterpolation(
        -1.0e6*rj.forcing_data[:, rj.pars.timecolumn.v],  Float64.(rj.forcing_data[:, rj.pars.datacolumn.v]), 
        extrapolation_bc = rj.pars.extrap_value.v )  # fill out of range values

    var_tforce = PB.VarDepScalar("tforce",     "yr",  
        "historical time at which to apply forcings, present = 0 yr")
    var_FORCE = PB.VarPropScalar(rj.pars.forcename.v,    "",  
        "forcing interpolated from $(forcingfile) column $(rj.pars.datacolumn.v)")

    PB.setfrozen!.(PB.get_parameters(rj)) # no modification after spreadsheet read

    PB.add_method_do!(
        rj, 
        do_force_spreadsheet,
        (PB.VarList_single(var_tforce), PB.VarList_single(var_FORCE)),
        p=rj.interp_FORCE
    )

    return nothing
end

function do_force_spreadsheet(m::PB.ReactionMethod, (var_tforce, var_FORCE), cellrange::PB.AbstractCellRange, deltat)
    interp_FORCE = m.p

    tforce = var_tforce[]
    var_FORCE[] =  interp_FORCE(tforce)
    
    return nothing
end



"Install create_reactionXXX factories when module imported"
function __init__()
    PB.add_reaction_factory(ReactionForce_CK_Solar)
    PB.add_reaction_factory(ReactionForce_UDWEbergman2004)
    PB.add_reaction_factory(ReactionForce_Bbergman2004)
    PB.add_reaction_factory(ReactionForce_CPlandrelbergman2004)
    PB.add_reaction_factory(ReactionForce_spreadsheet)   
    return nothing
end


end # module
