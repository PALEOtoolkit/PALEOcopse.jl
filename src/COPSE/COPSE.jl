# COPSE/COPSE.jl

module COPSE

function srcdir()
    return dirname(pathof(COPSE))
end


include("ModelBergman2004.jl")
include("MapAtmOceanReservoirs.jl")

include("CIsotopes.jl")


"""
    copse_crash(value, label, tmodel) -> outputvalue

COPSE 5_14 (Bergman 2004). Enforce +ve reservoir values by limiting fluxes at small reservoir sizes (??)
`output` = 1 ie independent of `value` for `value` > 0.1 (implicitly the usual range), 
but limited towards zero for `value` < 0.1
"""
function copse_crash(value, label, tmodel)
    
    valuelowlimit = 0.1
    if (value > valuelowlimit)
        output = 1.0
    elseif (value > 0.0)
        output = value*1.0/valuelowlimit
        @warn "copse_crash $label value $value tmodel $tmodel"
    else
        output = 0.0
        @warn "copse_crash $label value $value tmodel $tmodel"
    end
    return output
end

end # modules
