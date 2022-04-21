

using Test
using Documenter


@testset "PALEOcopse all" begin

@testset "PALEOcopse" begin
include("forcings/runforcingtests.jl")

include("global/runtemperaturetests.jl")

doctest(PALEOcopse; manual=false)  

end 

@testset "PALEOcopse/examples/COPSE" begin

include("../examples/COPSE/runtests.jl")

end

end