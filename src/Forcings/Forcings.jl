
module Forcings

function srcdir()
    return @__DIR__
end


include("COPSEForcings.jl")

include("LipForcing.jl")

end