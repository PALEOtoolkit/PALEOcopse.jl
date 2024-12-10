using Logging
import DataFrames

using Plots

import PALEOboxes as PB
import PALEOmodel

import PALEOcopse
import PALEOocean

include("copse_reloaded_bergman2004_expts.jl")
include("compare_output.jl")


global_logger(ConsoleLogger(stderr,Logging.Info))

# logfile = open("COPSE_reloaded_bergman2004_log.txt", "w")
# println("writing log output to ", logfile)
# global_logger(SimpleLogger(logfile, Logging.Info))


@info "Start $(@__FILE__)"



# Baseline Phanerozoic configuration
# comparemodel = CompareOutput.copse_output_load("bergman2004","")
# comparemodel = CompareOutput.copse_output_load("reloaded","original_baseline")
comparemodel = nothing
model = PB.create_model_from_config(
    joinpath(@__DIR__, "COPSE_reloaded_bergman2004_cfg.yaml"), "model1"
)
copse_reloaded_bergman2004_expts(model, ["baseline"])
tspan=(-1000e6, 0)


# No S cycle
# comparemodel = CompareOutput.copse_output_load("reloaded","original_baseline")
# model = PB.create_model_from_config(
#     joinpath(@__DIR__, "COPSE_reloaded_bergman2004_cfg.yaml"), "model1";
#     modelpars=Dict("enableS"=>false, "disableS"=>true),
# )
# 
# tspan=(-1000e6, 0)


# Modern steady-state (constant forcings) with CO2 pulse
# comparemodel=nothing
# model = PB.create_model_from_config(
#     joinpath(@__DIR__, "COPSE_reloaded_bergman2004_cfg.yaml"), "model1";
#     modelpars=Dict("tforcevar"=>"global.tforce_constant"),
# )
# copse_reloaded_bergman2004_expts(
#     model, 
#     [
#         ("tforce_constant", 0.0),
#         "VCI",
#         "mocbProdLinear",
#         "noNcycle",
#         ("CO2pulse", 6e18, 0.1e6, 1e5),
#         ("setpar", "global", "reservoir_A", "f_atfrac", "quadratic"),
#     ], 
# )
# tspan=(-10e6, 10e6)

initial_state, modeldata = PALEOmodel.initialize!(model)

# call ODE function to check derivative
initial_deriv = similar(initial_state)
PALEOmodel.SolverFunctions.ModelODE(modeldata)(initial_deriv, initial_state , nothing, 0.0)
println("initial_state", initial_state)
println("initial_deriv", initial_deriv)

println("integrate, no jacobian")
                                          
paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

# NB first run includes JIT time
@time PALEOmodel.ODE.integrate(
    paleorun, initial_state, modeldata, tspan,
    solvekwargs=( # https://diffeq.sciml.ai/stable/basics/common_solver_opts/
        reltol=1e-4,
        # saveat=1e6,
        dtmax=1e5,  # needed for short pulses
    )
)



##############################
# Plot output
###############################

# no subplots
# pager=PALEOmodel.DefaultPlotPager()

# assemble subplots onto screens
# plotlyjs(size=(1200, 700)) # useless - PlotlyJS backend can't handle per-subplot legends !
gr(size=(1600, 900)) # GR backend with large screen size
pager=PALEOmodel.PlotPager(9, (legend_background_color=nothing, ))


copse_reloaded_bergman2004_plot(paleorun.output, pager=pager)

###################################
# Check diff against archived Matlab model output
######################################

if !isnothing(comparemodel)
    diff = CompareOutput.compare_copse_output(paleorun.output, comparemodel)
    firstpoint=100
    show(sort(DataFrames.describe(diff[firstpoint:end, :], :std, :min, :max, :mean), :std), allrows=true)

    println()

    # Plots to compare two model runs
    CompareOutput.copse_reloaded_comparecopse(paleorun.output, comparemodel, include_Sr=false, mccb_extrafields=true, pager=pager) 
end

####################################################################################################
# Overlay comparison data (data compilations available on request, not publically distributable)
#################################################################################################

if false
    # include("compare_data.jl")
    # copse_reloaded_comparedata([paleorun.output], include_Sr=false, pager=pager)
else
    copse_reloaded_bergman2004_plot_summary([paleorun.output], pager=pager)
end


@info "End $(@__FILE__)"

if @isdefined logfile
    close(logfile)
end
