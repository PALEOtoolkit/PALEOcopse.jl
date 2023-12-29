using Test
using NBInclude
using Logging
using Plots
import DataFrames

import PALEOboxes as PB
import PALEOmodel

import PALEOcopse

include("compare_output.jl")

@testset "COPSE examples" begin
skipped_testsets = [
    # "COPSE_reloaded_reloaded",
    # "COPSE_reloaded_bergman2004",
    # "COPSE_bergman2004_bergman2004",
    # "COPSE_reloaded_reloaded.ipynb",
]

!("COPSE_reloaded_reloaded" in skipped_testsets) &&
@testset "COPSE_reloaded_reloaded" begin

    include("copse_reloaded_reloaded_expts.jl")

    # load archived model output
    comparemodel = CompareOutput.copse_output_load("reloaded", "reloaded_baseline")

    model = copse_reloaded_reloaded_expts("reloaded", ["baseline"])
    initial_state, modeldata = PALEOmodel.initialize!(model)

    run = PALEOmodel.Run(model=model, output=PALEOmodel.OutputWriters.OutputMemory())

    @time PALEOmodel.ODE.integrate(
        run, initial_state, modeldata, (-1000e6, 0),
        solvekwargs=(
            dtmin=0.0,
            reltol=1e-4,
        )
    )
    # @time PALEOmodel.ODE.integrateForwardDiff(run, initial_state, modeldata, (-1000e6, 0), jac_ad_t_sparsity=0.0, solvekwargs=(reltol=1e-4,))

    # conservation checks
    conschecks = [
        ("total_C",     :v,             1e-6 ),
        ("total_C",     :v_moldelta,    1e-6),
        ("total_S",     :v,             1e-6),
        ("total_S",     :v_moldelta,    1e-6),
        ("total_redox", nothing,        1e-6)
    ]
    for (varname, propertyname, rtol) in conschecks
        startval, endval = PB.get_property(
            PB.get_data(run.output, "global."*varname),
            propertyname=propertyname,
        )[[1, end]]
        println("check $varname $startval $endval $rtol")
        @test isapprox(startval, endval, rtol=rtol)
    end

    # comparison to archived output
    diff = CompareOutput.compare_copse_output(run.output, comparemodel)
    firstpoint=50
    diffsummary = DataFrames.describe(diff[firstpoint:end, :], :std, :min, :max, :mean)

    stdlimits = [
        (:mccb_delta_diff,      1.5e-3),
        (:Sr_delta_diff,        4.0e-7),
        (:pCO2PAL_reldiff,      6.0e-4),
        (:pO2PAL_reldiff,       3.0e-4),
        (:S_delta_diff,         2.5e-3),
        (:P_reldiff,            5.0e-4),
        (:ANOX_diff,            3.0e-3),
    ]
    for (var, limit) in stdlimits
        print(var)
        vardiff = diffsummary[diffsummary.variable.==var, "std"][]
        println("  diff=", vardiff, "  limit=", limit)
        @test diffsummary[diffsummary.variable.==var, "std"][] < limit
    end

end

!("COPSE_reloaded_bergman2004" in skipped_testsets) &&
@testset "COPSE_reloaded_bergman2004" begin
    include("copse_reloaded_bergman2004_expts.jl")

    comparemodel = CompareOutput.copse_output_load("reloaded","original_baseline")

    model = copse_reloaded_bergman2004_expts("bergman2004", ["baseline"])

    initial_state, modeldata = PALEOmodel.initialize!(model)

    run = PALEOmodel.Run(model=model, output=PALEOmodel.OutputWriters.OutputMemory())

    @time PALEOmodel.ODE.integrate(
        run, initial_state, modeldata, (-1000e6, 0),
        solvekwargs=(
            reltol=1e-4,
            dtmin=0.0,
            # saveat=1e6,
        ),
    )

    # conservation checks
    conschecks = [
        ("total_C",     :v,             1e-6 ),
        ("total_C",     :v_moldelta,    1e-6),
        ("total_S",     :v,             1e-6),
        ("total_S",     :v_moldelta,    1e-6),
        ("total_redox", nothing,        1e-6)
    ]
    for (varname, propertyname, rtol) in conschecks
        startval, endval = PB.get_property(
            PB.get_data(run.output, "global."*varname),
            propertyname=propertyname
        )[[1, end]]
        println("check $varname $startval $endval $rtol")
        @test isapprox(startval, endval, rtol=rtol)
    end

    # comparison to archived output
    diff = CompareOutput.compare_copse_output(run.output, comparemodel)
    firstpoint = 100
    diffsummary = DataFrames.describe(diff[firstpoint:end, :], :std, :min, :max, :mean)

    stdlimits = [
        (:mccb_delta_diff,      4.0e-3),
        (:pCO2PAL_reldiff,      3.0e-4),
        (:pO2PAL_reldiff,       4.0e-4),
        (:S_delta_diff,         6.0e-3),
        (:P_reldiff,            2.0e-4),
        (:ANOX_diff,            2.0e-4),
    ]
    for (var, limit) in stdlimits
        print(var)
        vardiff = diffsummary[diffsummary.variable.==var, "std"][]
        println("  diff=", vardiff, "  limit=", limit)
        @test diffsummary[diffsummary.variable.==var, "std"][] < limit
    end

end

!("COPSE_bergman2004_bergman2004" in skipped_testsets) &&
@testset "COPSE_bergman2004_bergman2004" begin
    include("copse_bergman2004_bergman2004_expts.jl")

    comparemodel = CompareOutput.copse_output_load("bergman2004","")

    model = copse_bergman2004_bergman2004_expts(["baseline"])

    initial_state, modeldata = PALEOmodel.initialize!(model)

    run = PALEOmodel.Run(model=model, output=PALEOmodel.OutputWriters.OutputMemory())

    @time PALEOmodel.ODE.integrate(
        run, initial_state, modeldata, (-600e6, 0),
        solvekwargs=(
            dtmin=0.0,
            saveat=1e6,
        )
    )

    # conservation checks
    conschecks = [
        ("total_C",     :v,             1e-6 ),
        ("total_C",     :v_moldelta,    1e-6),
        ("total_S",     :v,             1e-6),
        ("total_S",     :v_moldelta,    1e-6),
        ("total_redox", nothing,        1e-6)
    ]
    for (varname, propertyname, rtol) in conschecks
        startval, endval = PB.get_property(
            PB.get_data(run.output, "global."*varname),
            propertyname=propertyname
        )[[1, end]]
        println("check $varname $startval $endval $rtol")
        @test isapprox(startval, endval, rtol=rtol)
    end

    # comparison to archived output
    diff = CompareOutput.compare_copse_output(run.output, comparemodel)
    firstpoint = 1
    diffsummary = DataFrames.describe(diff[firstpoint:end, :], :std, :min, :max, :mean)

    stdlimits = [
        (:mccb_delta_diff,      3.0e-2),
        (:pCO2PAL_reldiff,      9.0e-3),
        (:pO2PAL_reldiff,       6.0e-3),
        (:S_delta_diff,         4.0e-2),
        (:P_reldiff,            6.0e-3),
        (:ANOX_diff,            2.0e-3),
    ]
    for (var, limit) in stdlimits
        print(var)
        vardiff = diffsummary[diffsummary.variable.==var, "std"][]
        println("  diff=", vardiff, "  limit=", limit)
        @test diffsummary[diffsummary.variable.==var, "std"][] < limit
    end

end

!("COPSE_reloaded_reloaded.ipynb" in skipped_testsets) &&
@testset "COPSE_reloaded_reloaded.ipynb" begin
    # unicodeplots()
    gr() # headless environment will require ENV["GKSwstype"] = "100"
    @nbinclude("COPSE_reloaded_reloaded.ipynb"; counters=2:1000) # omit first cell with display etc setup
    print("\n")
end

end
