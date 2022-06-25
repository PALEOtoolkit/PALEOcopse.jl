module LipForcing

import Interpolations

import XLSX  # Excel file access
import Roots

import PALEOboxes as PB
using PALEOboxes.DocStrings

import PALEOcopse

import Infiltrator # Julia debugger

##################################################
# Forcing due to one LIP
##################################################

"forcing due to one LIP"
Base.@kwdef struct ForceLIP
    Name::String
    Type::String
    "emplacement time, yr relative to present day"
    liptime::Float64
    "km^2 peak area exposed to subaerial weathering"
    peakCFBarea::Float64
    "km^2 present-day area, or missing if not available"
    presentCFBarea::Union{Missing,Float64}
    "co2 released (mol C)"
    co2_potential::Float64
    "co2 release isotopic composition"
    co2_d13C::Float64   = -5.0
    "yr wait before erosion begins"
    decayoffset::Float64
    "1/yr lambda for exponential area decay"
    lamb::Float64
    "yr timeframe for CO2 release (full width of gaussian)"
    co2releasetime::Float64
    "true for sigmoidal emplacement"
    smoothempl::Bool
end

"create a sigmoidal or step rise and exponential decline in area"
function calc_LIP_CFBarea(lip::ForceLIP, timeyr)
    function sigmoid(smooth, t)
        s = 1.0/(1+exp(-t*smooth))
        return s
    end

    smoothing = 3e-6 # smooth curve to build up area
    offset = 1.5e6 # shift curve
            
    if lip.smoothempl
        # emplacement function is made from sigmoid
        empl = sigmoid(smoothing, timeyr - lip.liptime + offset)
    else
        # step at obj.liptime
        empl = 0.5 + 0.5.*sign(timeyr - lip.liptime)
    end
    # exponential decay
    erode =  (( 0.5 + 0.5*sign(timeyr - lip.liptime - lip.decayoffset) )
                        * ( 1 - exp(-lip.lamb*(timeyr - lip.liptime - lip.decayoffset) ) ))
            
    area =  lip.peakCFBarea*(empl - erode)
end

"calculate CO2 release"
function calc_LIP_co2(lip::ForceLIP, CIsotopeType, timeyr)
    function gaussnorm(width,offset)
        g = 1.0/(width*(2*pi)^0.5)*exp(-offset^2/(2*width^2))
        return g
    end
    
    co2release =  lip.co2_potential*gaussnorm(0.5*lip.co2releasetime, timeyr - lip.liptime) 

    return @PB.isotope_totaldelta(CIsotopeType, co2release, lip.co2_d13C)
end


######################################################
# Create set of LIPs from definitions in spreadsheet
#####################################################

"""
    read_lips_xlsx(lipfile)  -> LIPdata::NamedTuple

read LIP data from spreadsheet
"""
function read_lips_xlsx(lipfile)

    @info "    read_lips_xlsx: spreadsheet $(lipfile)"

    xf = XLSX.readxlsx(lipfile)
    sh = xf["Sheet1"]

    column_names = sh["A1:I1"]
    @info "    read_lips_xlsx: column_names: $(column_names)"
    expected_column_names = ["Age" "Name" "Type" "Continental areas" "All volumes" "CO2 release (min)" "CO2 release (max)" "Degassing duration" "Present day area"]
    column_names == expected_column_names || error("spreadsheet column names don't match expected_column_names=", expected_column_names)

    Age             = Float64[]
    Name            = String[]
    Type            = String[]
    CFBarea         = Union{Missing, Float64}[]
    Allvolume       = Float64[]
    CO2min          = Float64[]
    CO2max          = Float64[]
    DegassDuration  = Float64[]
    PresentDayArea  = Union{Missing,Float64}[]

    current_row = 3  # skip two header rows
    while !ismissing(sh[current_row, 1]) # read rows until first column (Age) has a missing value
        push!(Age,              sh[current_row, 1])
        push!(Name,             sh[current_row, 2])
        push!(Type,             sh[current_row, 3])
        cfb_area              = sh[current_row, 4]
        push!(CFBarea,      cfb_area == "-" ? missing : cfb_area)
        push!(Allvolume,        sh[current_row, 5])
        push!(CO2min,           sh[current_row, 6])
        push!(CO2max,           sh[current_row, 7])
        push!(DegassDuration,   sh[current_row, 8])
        # push!(PresentDayArea,   sh[current_row, 9])
        pda = sh[current_row,9]        
        if pda isa String
            @warn "read_lips_xlsx: row $current_row unexpected String value for 'Present day area' $pda (replacing with missing)"
            pda = missing
        end
        push!(PresentDayArea,   pda)
        current_row += 1
    end

    return (;Age, Name, Type, CFBarea, Allvolume, CO2min, CO2max, DegassDuration, PresentDayArea) # NamedTuple
end

"""
    create_LIPs(LIPdata, co2releasefield::Symbol, default_lambda; [lipkwargs...]) -> LIPs::Vector{ForceLIP}

create set of LIPs from spreadsheet data
"""
function create_LIPs(
    LIPdata, co2releasefield::Symbol, default_lambda; 
    LIPtrange=(-Inf,Inf),
    smoothempl=false,
    areafac=1.0,
    decayoffset=0.0,
    default_co2releasetime=2e5
)

    LIPs = ForceLIP[]

    function missingtodefault(val, defaultval)
        if ismissing(val)
            return defaultval
        else
            return val
        end
    end

    for i in 1:length(LIPdata.Age)
        liptime = -1.0e6*LIPdata.Age[i]
        if liptime >= LIPtrange[1] && liptime <= LIPtrange[2]
            if co2releasefield==:NoCO2
                lipCO2 = 0.0
            else
                lipCO2 = getfield(LIPdata, co2releasefield)[i]
            end
            
            co2releasetime = missingtodefault(1e6*LIPdata.DegassDuration[i], default_co2releasetime)
            
            peak_area   = areafac*missingtodefault(LIPdata.CFBarea[i], 0.0)

            present_area = LIPdata.PresentDayArea[i]
            if ismissing(present_area)
                lambda = default_lambda
            else
                # calculate decay rate to give present-day areafac
                lambda = log(present_area/peak_area)/(liptime+decayoffset)
            end

            push!(
                LIPs, 
                ForceLIP(
                    Name=LIPdata.Name[i],
                    Type=LIPdata.Type[i],
                    liptime=liptime,
                    peakCFBarea=peak_area,
                    presentCFBarea=present_area,
                    co2_potential=lipCO2,
                    decayoffset=decayoffset,
                    lamb=lambda,
                    co2releasetime=co2releasetime,
                    smoothempl=smoothempl
                )
            )
        end
    end
        
    return LIPs
end


"solve for default_lambda required for set of LIPs to give present_CFB_area_target"
function find_default_lambda(present_CFB_area_target, LIPdata, lipkwargs )
    "calculate present_CFB_area given default_lambda"
    function present_CFB_area_error(default_lambda)
        LIPs = create_LIPs(LIPdata, :NoCO2, default_lambda; lipkwargs...)
        pCFBa = 0.0
        for lip in LIPs
            pCFBa += calc_LIP_CFBarea(lip, 0.0)
        end
        return pCFBa - present_CFB_area_target
    end

    return Roots.find_zero(present_CFB_area_error, (0.25e-8, 4e-8), Roots.Bisection())
end

"""
    ReactionForce_LIPs

LIP basalt area and optional CO2 forcing ([Mills2014b](@cite)), updated in [Lenton2018](@cite) COPSE Reloaded).

Reads spreadsheet table of LIP data, and creates per-LIP area vs time profile including weathering reduction in area.
A default area reduction rate (applied to those LIPs where it is poorly constrained) is tuned to match overall modern continental
flood basalt area.

If Parameter `co2releasefield` != `NoCO2`, a CO2 pulse (duration 2e5yr, magnitude from spreadsheet, -5 per mil) is output as
`LIP_CO2` and added to Variable `flux_C` -> `fluxSedCrusttoAOcean.flux_C`.

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionForce_LIPs{P} <:  PB.AbstractReaction
    base::PB.ReactionBase
    
    pars::P = PB.ParametersTuple(
        PB.ParString("datafolder", PALEOcopse.Forcings.srcdir(),
            description="folder for forcing data 'datafile'"),
        PB.ParString("datafile", "CR_force_LIP_table.xlsx", 
            description="spreadsheet with forcing data (path relative to 'datafolder')"),    
        PB.ParDouble("present_day_CFB_area", 4.8e6, units="km^2", 
            description="present day continental flood basalt area"),
        PB.ParBool("smoothempl", false, 
            description="smooth LIP basalt emplacement"),
        PB.ParString("co2releasefield", "NoCO2", allowed_values=["NoCO2", "CO2min", "CO2max"],
            description="LIP CO2 release"),
        PB.ParType(PB.AbstractData, "CIsotope", PB.ScalarData,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable carbon isotopes and specify isotope type"),
    )

    LIP_data = nothing
    LIPs::Vector{ForceLIP} = ForceLIP[]

    default_lambda  = nothing

end


function PB.register_methods!(rj::ReactionForce_LIPs)
     
    CIsotopeType = rj.pars.CIsotope.v
    PB.setfrozen!(rj.pars.CIsotope)
    @info "register_methods! ReactionForce_LIPs: $(PB.fullname(rj)) CIsotopeType=$(CIsotopeType)"

    vars = [
        PB.VarDepScalar("tforce",     "yr",  "historical time at which to apply forcings, present = 0 yr"),
        PB.VarPropScalar("CFB_area",  "km^2",  "LIP continental flood basalt area"),
        PB.VarPropScalar("LIP_CO2",  "mol yr-1",  "LIP CO2 release",
            attributes=(:field_data=>CIsotopeType,)),
        PB.VarContribScalar("flux_C"=>"fluxSedCrusttoAOcean.flux_C",  "mol yr-1",  "LIP CO2 release", 
            attributes=(:field_data=>CIsotopeType,)),
    ]
    
    PB.add_method_setup!(
        rj,
        setup_force_LIPs,
        (),
    )

    PB.add_method_do!(
        rj,
        do_force_LIPs,
        (PB.VarList_namedtuple(vars), ), 
        p=CIsotopeType,
    )

    return nothing
end

function setup_force_LIPs(m::PB.ReactionMethod, (), cellrange::PB.AbstractCellRange, attribute_name)
    rj = m.reaction

    attribute_name == :setup || return nothing

    lipfile = joinpath(rj.pars.datafolder.v, rj.pars.datafile.v)
    @info "setup_force_LIPs! ReactionForce_LIPs: $(PB.fullname(rj)) loading LIP forcing from 'datafolder/datafile'='$(lipfile)'"

    rj.LIP_data = read_lips_xlsx(lipfile)

    lipkwargs = (smoothempl=rj.pars.smoothempl.v, )

    rj.default_lambda = find_default_lambda(rj.pars.present_day_CFB_area.v, rj.LIP_data, lipkwargs )

    @info "    found default_lambda=$(rj.default_lambda) "*
        "to match present_day_CFB_area = $(rj.pars.present_day_CFB_area.v) km^2"
 
    co2releasefield = Symbol(rj.pars.co2releasefield.v)
    @info "    add LIP CO2 release from $(rj.pars.co2releasefield.v) field"

    empty!(rj.LIPs)
    rj.LIPs = create_LIPs(rj.LIP_data, co2releasefield, rj.default_lambda; lipkwargs...)
        
    return nothing
end


function do_force_LIPs(m::PB.ReactionMethod, (vars, ), cellrange::PB.AbstractCellRange, deltat)
    rj = m.reaction
    CIsotopeType = m.p

    vars.CFB_area[] = 0.0
    for lip in rj.LIPs
        vars.CFB_area[] += calc_LIP_CFBarea(lip, PB.value_ad(vars.tforce[]))
    end

    vars.LIP_CO2[] = @PB.isotope_totaldelta(CIsotopeType, 0.0, 0.0)    
    for lip in rj.LIPs
        lip_CO2 = calc_LIP_co2(lip, CIsotopeType, PB.value_ad(vars.tforce[]))
        vars.LIP_CO2[] += lip_CO2
        vars.flux_C[] += lip_CO2
    end
   
    return nothing
end



end # module
