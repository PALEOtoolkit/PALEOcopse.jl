

# COPSE Reloaded

These examples demonstrate COPSE reloaded ([Lenton2018](@cite)), modularised with an atm/ocean/land/sedcrust split.

## Activate the PALEOcopse/examples environment

First, change the Julia REPL working directory to the `PALEOcopse/examples/COPSE` folder:

In `VS code`, right click on this folder in the file browser and select `Julia: Change to This Directory`. Or from the REPL, use the `cd` command):

    julia> cd("examples/COPSE")

If it is not already active, activate the Julia environment `PALEOcopse/examples`:

In `VS code`, right click on `PALEOjulia/PALEOexamples` or any subfolder in the file browser and select `Julia: Activate Parent Environment`. Or from the REPL, use `]` to enter package management:

```julia
pwd()
"/home/sd336/software/julia/PALEOcopse/examples/COPSE"
import Pkg
Pkg.activate("../")
```

## To run the COPSE Reloaded baseline example for the Phanerozoic:
   
    julia> include("COPSE_reloaded_reloaded.jl")

This will run and plot output (NB: the first run will be slow as Julia JIT compiles the code).

## To display model Parameters, Variables, and output.

See [PALEOmodel](https://github.com/PALEOtoolkit/PALEOmodel.jl) [documentation](https://paleotoolkit.github.io/PALEOmodel.jl)

## Change climate sensitivity from default 3C to 6C and compare pCO2 prediction
See [Lenton2018](@cite) Fig. 11

Modify parameters and rerun model:
```julia
default_3C_output = paleorun.output  # keep default 3C output
PB.setvalue!(PB.get_reaction(model, "global", "temp_global").pars.k_c, 8.656)  # change climate sensitivity
paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())
PALEOmodel.ODE.integrate(paleorun, initial_state, modeldata, (-1000e6, 0), solvekwargs=(reltol=1e-4,)); # rerun
```

Compare pCO2:
```julia
using Plots; plotlyjs(size=(750, 500))  # load Julia Plots.jl and choose PlotlyJS backend
pCO2atm_3C = PALEOmodel.get_array(default_3C_output, "atm.pCO2atm") # get FieldArray for plotting
plot(1e6*pCO2atm_3C, label="3C default", ylabel="pCO2 (1e-6 atm)") # plot x1e6
pCO2atm_6C = PALEOmodel.get_array(paleorun.output, "atm.pCO2atm")   # get FieldArray
plot!(1e6*pCO2atm_6C, label="6C", ylabel="pCO2 (1e-6 atm)") # plot x1e6
plot!(xlim=(-600e6, 0), ylim=(0, 4000))  # rescale axes
```

## Jupyter notebook version
The same content as above (with additional examples) is available as the Jupyter notebook `COPSE_reloaded_reloaded.ipynb`.

To start Jupyter from within the Julia REPL:

    julia> using IJulia
    julia> notebook(dir=pwd(), detached=true)