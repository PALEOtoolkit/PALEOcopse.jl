using Logging
import DataFrames

using Plots 

import PALEOboxes as PB
import PALEOmodel

import PALEOcopse

include("copse_reloaded_reloaded_expts.jl")
include("compare_output.jl")

global_logger(ConsoleLogger(stderr,Logging.Info))

# logfile = open("COPSE_reloaded_reloaded_log.txt", "w")
# println("writing log output to ", logfile)
# global_logger(SimpleLogger(logfile, Logging.Info))

@info "Start $(@__FILE__)"

# load archived model output 
# comparemodel=nothing


# Baseline Phanerozoic configuration
# comparemodel = CompareOutput.copse_output_load("reloaded","reloaded_baseline")
comparemodel = nothing
run = copse_reloaded_reloaded_expts(
    "reloaded",
    ["baseline"],
)
tspan=(-1000e6, 0)


# OOE oscillations with VCI and linear mocb
# comparemodel = CompareOutput.copse_output_load("reloaded","reloaded_baseline")
# run = copse_reloaded_reloaded_expts(
#     "reloaded", 
#     [
#         "VCI",
#         "mocbProdLinear",
#         ("setpar", "global", "force_LIPs", "co2releasefield", "CO2max"),
#     ],
# ) 
# tspan=(-1000e6, 1e6)

######################################
# Integrate model
###########################################

initial_state, modeldata = PALEOmodel.initialize!(run)

# call ODE function to check derivative
initial_deriv = similar(initial_state)
PALEOmodel.ODE.ModelODE(modeldata)(initial_deriv, initial_state , (run=run, modeldata=modeldata), 0.0)
println("initial_state", initial_state)
println("initial_deriv", initial_deriv)

println("integrate, no jacobian")
@time PALEOmodel.ODE.integrate(
    run, initial_state, modeldata, tspan, 
    solvekwargs=( # https://diffeq.sciml.ai/stable/basics/common_solver_opts/
        reltol=1e-4,
        # reltol=1e-5,
        # saveat=1e6,
        # dtmax=1e4,
    )
) 

# println("integrate, jacobian from autodifferentiation")
# @time PALEOmodel.ODE.integrateForwardDiff(run, initial_state, modeldata, (-1000e6, 0), jac_ad_t_sparsity=0.0, solvekwargs=(reltol=1e-5,))
                                   
# PB.TestUtils.bench_model(run.model)



##############################
# Plot output
###############################

# no subplots
# pager=PALEOmodel.DefaultPlotPager()

# assemble subplots onto screens
# plotlyjs(size=(1200, 700)) # useless - PlotlyJS backend can't handle per-subplot legends !
gr(size=(1600, 900)) # GR backend with large screen size
pager=PALEOmodel.PlotPager(9, (legend_background_color=nothing, ))


copse_reloaded_reloaded_plot(run.output, pager=pager)

#####################################
# Check diff against comparison model 
########################################

if !isnothing(comparemodel)   
    diff = CompareOutput.compare_copse_output(run.output, comparemodel)
    firstpoint=50
    show(sort(DataFrames.describe(diff[firstpoint:end, :], :std, :min, :max, :mean), :std), allrows=true)

    println()

    # Overlay plots for comparison
    CompareOutput.copse_reloaded_comparecopse(run.output, comparemodel, include_Sr=true, pager=pager) 
end

####################################################################################################
# Overlay comparison data (data compilations available on request, not publically distributable)
#################################################################################################

if false
    include("compare_output_plots.jl")
    copse_reloaded_comparedata([run.output], include_Sr=true, pager=pager)
else
    copse_reloaded_reloaded_plot_summary([run.output], pager=pager)
end

pager(:newpage) # flush output

@info "End $(@__FILE__)"

if @isdefined logfile
    close(logfile)
end
