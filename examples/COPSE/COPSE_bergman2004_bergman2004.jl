using Logging
using DataFrames
using Plots

global_logger(ConsoleLogger(stderr,Logging.Info))

@info "Start $(@__FILE__)"

@info "importing modules ... (may take a few seconds)"
import PALEOboxes as PB
import PALEOmodel

import PALEOcopse
@info "                  ... done"

include("compare_output.jl")
include("copse_bergman2004_bergman2004_expts.jl")


# load archived model output
comparemodel = CompareOutput.copse_output_load("bergman2004","")
# comparemodel = nothing

use_TEMP_DAE = false

model = copse_bergman2004_bergman2004_expts(
    ["baseline"],
    modelpars=Dict("temp_DAE"=>use_TEMP_DAE),
)

initial_state, modeldata = PALEOmodel.initialize!(model)

# call ODE function to check derivative
initial_deriv = similar(initial_state)
PALEOmodel.SolverFunctions.ModelODE(modeldata)(initial_deriv, initial_state , nothing, 0.0)
println("initial_state", initial_state)
println("initial_deriv", initial_deriv)

run = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

if use_TEMP_DAE
    println("integrate, DAE")
    # first run includes JIT time
    @time PALEOmodel.ODE.integrateDAE(
        run, initial_state, modeldata, (-600e6, 0),
        solvekwargs=(
            dtmin=0.0,
            saveat=1e6, # save output every 1e6 yr see https://diffeq.sciml.ai/dev/basics/common_solver_opts/
        )
    )

    # sparse jacobian with IDA(linear_solver=:KLU) fails at ~-350e6 yr ?
    # @time PALEOmodel.ODE.integrateDAEForwardDiff(run, initial_state, modeldata, (-650e6, 0), alg=IDA(linear_solver=:KLU))
else
    println("integrate, ODE")
    # first run includes JIT time
    @time PALEOmodel.ODE.integrate(
        run, initial_state, modeldata, (-600e6, 0),
        solvekwargs=(
            dtmin=0.0,
            saveat=1e6, # save output every 1e6 yr see https://diffeq.sciml.ai/dev/basics/common_solver_opts/
        )
    )
    # sparse jacobian with IDA(linear_solver=:KLU) fails at ~-350e6 yr ?
    # @time PALEOmodel.ODE.integrateForwardDiff(run, initial_state, modeldata, (-650e6, 0), alg=CVODE_BDF(linear_solver=:KLU))
end


if !isnothing(comparemodel)
    diff = CompareOutput.compare_copse_output(run.output, comparemodel)
    firstpoint=1
    show(sort(DataFrames.describe(diff[firstpoint:end, :], :std, :min, :max, :mean), :std), allrows=true)
end
println()


##############################
# Plot output
###############################

# no subplots
# pager=PALEOmodel.DefaultPlotPager()

# assemble subplots onto screens
gr(size=(1200, 700)) # GR backend with large screen size
pager=PALEOmodel.PlotPager(6, (legend_background_color=nothing, ))

copse_bergman2004_bergman2004_plot(run.output; pager=pager)
if !isnothing(comparemodel)
    copse_bergman2004_bergman2004_compare([run.output, comparemodel]; pager=pager)
end


@info "End $(@__FILE__)"
