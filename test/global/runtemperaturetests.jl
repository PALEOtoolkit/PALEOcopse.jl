
using Test

import PALEOboxes as PB
import PALEOmodel
import PALEOcopse

@testset "GlobalTemperature" begin

    model = PB.create_model_from_config(joinpath(@__DIR__, "configtemperature.yaml"), "model1")

    # Test domains
    @test PB.get_num_domains(model) == 1
    global_domain = PB.get_domain(model,"global")   
    @test PB.get_length(global_domain) == 1

    modeldata = PB.create_modeldata(model)
    PB.allocate_variables!(model, modeldata, hostdep=false)
    @test length(PB.get_unallocated_variables(global_domain, modeldata)) == 3
    @test PB.check_ready(model, modeldata, throw_on_error=false) == false
    # allocate state and sms variables
    PB.allocate_variables!(model, modeldata, hostdep=true)   
    PALEOmodel.set_default_solver_view!(model, modeldata)
    @test PB.check_ready(model, modeldata) == true

    # get modelcreated variables
    modelcreated_data = Dict(var.name=>PB.get_data(var, modeldata) for var in PB.get_variables(global_domain, hostdep=false))
    # get host dependent variables
    hostdep_vars_vec = PB.get_variables(global_domain, hostdep=true)
    @test length(hostdep_vars_vec) == 3
    hostdep_data = Dict(var.name=>PB.get_data(var, modeldata) for var in hostdep_vars_vec)

    PB.initialize_reactiondata!(model, modeldata)
 
    # check Reaction configuration
    PB.check_configuration(model)
    # Initialise Reactions and non-state variables
    PB.dispatch_setup(model, :setup, modeldata)
    # Initialise state variables to norm_value
    PB.dispatch_setup(model, :norm_value, modeldata)
    PALEOmodel.copy_norm!(modeldata.solver_view_all)
    # Initialise state variables etc     
    PB.dispatch_setup(model, :initial_value, modeldata)

    println("after dispatch_setup - host-dependent variables:\n", hostdep_data)   

    # set to present-day
    hostdep_data["tforce"][] = 0.0

    PB.do_deriv(modeldata.dispatchlists_all)
    println("after do_deriv - model-created variables:\n", modelcreated_data)
    println("after do_deriv - host-dependent variables:\n", hostdep_data)

    # check that temperature correction is small
    @test abs(hostdep_data["TEMP_sms"][]) < 1.5e-2

    model = nothing
end
