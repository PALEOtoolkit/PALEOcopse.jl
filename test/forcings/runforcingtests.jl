
using Test
using NBInclude
using Plots

import PALEOboxes as PB

import PALEOcopse

@testset "Forcings" begin

@testset "Forcings base" begin

    model = PB.create_model_from_config(joinpath(@__DIR__, "configforcings.yaml"), "model1")

    # Test domains
    @test PB.get_num_domains(model) == 2

    global_domain = PB.get_domain(model, "global")
    @test global_domain.name == "global"
    @test PB.get_length(global_domain) == 1

    modeldata = PB.create_modeldata(model)
    PB.allocate_variables!(model, modeldata, hostdep=false)
    @test PB.check_ready(model, modeldata, throw_on_error=false) == false

    # state and sms variables
    PB.allocate_variables!(model, modeldata, hostdep=true)
    @test PB.check_ready(model, modeldata) == true
    PB.set_default_solver_view!(model, modeldata) # also (re)allocates tforce

    modelcreated_vars_dict = Dict([(var.name, var) for var in PB.get_variables(global_domain, hostdep=false)])

    PB.initialize_reactiondata!(model, modeldata)    
      
    @info "dispatch_setup"
    PB.dispatch_setup(model, :initial_value, modeldata)
   
    dispatchlists = modeldata.dispatchlists_all
 
    @info "do_deriv"

    if false # verbose
    for tforce in [-1e9 -250e6 0.0 1e9]
        PB.set_tforce!(modeldata.solver_view_all, tforce)
        PB.do_deriv(dispatchlists)
        println("tforce = ", tforce, " yr")
        for (name, var ) in modelcreated_vars_dict
            println("\t", PB.fullname(var), " = ", PB.get_data(var, modeldata))
        end
    end
    end

    PB.set_tforce!(modeldata.solver_view_all, 0.0)
    PB.do_deriv(dispatchlists)

    @test PB.get_data(modelcreated_vars_dict["SOLAR"], modeldata)[] == 1368.0
    @test PB.get_data(modelcreated_vars_dict["F_CONST"], modeldata)[] == 42.0

    model = nothing
end

if true
@testset "PlotForcings.ipynb" begin
    gr()  # headless environment will require ENV["GKSwstype"] = "100"
    @nbinclude("PlotForcings.ipynb"; counters=2:1000) # omit first cell with display etc setup
    println("done PlotForcings.ipynb")
end
end

end
