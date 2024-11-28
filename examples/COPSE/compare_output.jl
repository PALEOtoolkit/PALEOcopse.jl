module CompareOutput

import DataFrames
import Interpolations
import MAT
import CSV
import Plots

import PALEOboxes as PB
import PALEOmodel


"""
    copse_output_load(copse_version, copse_configuration) -> test_output::Dict{String, DataFrame}

Load archived Matlab COPSE model output, and convert format to Dict of DataFrames.

Matlab output is in test_output["global"].

Additional derived fields are in test_output["extrafields"].
"""
function copse_output_load(copse_version, copse_configuration)

    # Core tests output, available as part of repository / release package
    COPSE_test_output_dir = joinpath(splitdir(@__FILE__)[1], "COPSE_test_output")
    
    # Extra tests, available as a separate optional download
    COPSE_test_output_dir_full = joinpath(splitdir(@__FILE__)[1], "TODO")
    
    if copse_version == "bergman2004"
        test_output = copse_output_load_bergman2004(COPSE_test_output_dir)

    elseif copse_version == "lenton2016paleozoic"
        fname = "Lenton2016$(copse_configuration)"
        output_dirs = (COPSE_test_output_dir, COPSE_test_output_dir_full)
        folderpath = find_copse_output(output_dirs, fname);
        test_output = paleo_run.loadoutput(fname, folderpath)           
        # add field for delta_mccb
        test_output.diag.delta_mccb = test_output.diag.delta_A + test_output.diag.d_mccb;            
    elseif copse_version == "reloaded"
        output_dirs = (joinpath(COPSE_test_output_dir, "reloaded"), joinpath(COPSE_test_output_dir_full, "reloaded"))
        folderpath = find_copse_output(output_dirs, copse_configuration)
        test_output = paleo_matlab_output_load(copse_configuration, folderpath)
    
       
    else
        error("unrecognized copse_version ", copse_version)
    end

    return test_output
end

"Check a list of folder paths for specified fname.mat"
function find_copse_output(dirs, fname)
    
    folderpath = ""
    for dir in dirs
        testpath = joinpath(dir, fname*".mat")
        if isfile(testpath)
            folderpath = dir
            break
        end
    end

    if isempty(folderpath)
        error("test output file $(fname) not found")
    end

    return folderpath
end

"""
    paleo_matlab_output_load(copse_configuration, folderpath) -> test_output::PALEOmodel.OutputWriters::OutputMemory

Load data from PALEO Matlab run
"""
function paleo_matlab_output_load(filerootname, folderpath)

    outputfile = joinpath(folderpath, filerootname*".mat")

    @info "paleo_matlab_output_load loading output from file $outputfile"

    outputdata = MAT.matread(outputfile)

    T = outputdata["rs"]["T"][1,:] # Array size(1,n) -> Vector
    diag = outputdata["rs"]["diag"] # Arrays size(n, 1)
    S = outputdata["rs"]["S"]       # Arrays size(n, 1)  

    # dump everything into global domain
    glb = DataFrames.DataFrame()
    glb.tmodel = T

    for (key, value) in S
        DataFrames.insertcols!(glb, Symbol(key) => value[:,1])
    end
    for (key, value) in diag
        DataFrames.insertcols!(glb, Symbol(key) => value[:,1])
    end

     # add field for mccb_delta
    # DataFrames.insertcols!(glb, :mccb_delta => glb.delta_A + glb.d_mccb)   
    test_output_global = PALEOmodel.OutputWriters.OutputMemoryDomain("global", glb)

    # add fields for d13C
    extrafields = DataFrames.DataFrame()
    extrafields.tmodel = T
    DataFrames.insertcols!(extrafields, Symbol(:delta_atmos) => (glb.delta_A + glb.d_atmos))
    DataFrames.insertcols!(extrafields, Symbol(:delta_ocean) => (glb.delta_A + glb.d_ocean))
    DataFrames.insertcols!(extrafields, Symbol(:delta_mocb) => (glb.delta_A + glb.d_mocb))
    DataFrames.insertcols!(extrafields, Symbol(:delta_locb) => (glb.delta_A + glb.d_locb))
    DataFrames.insertcols!(extrafields, Symbol(:delta_mccb) => (glb.delta_A + glb.d_mccb))
    test_output_extrafields = 
        PALEOmodel.OutputWriters.OutputMemoryDomain("extrafields", extrafields)

    return PALEOmodel.OutputWriters.OutputMemory(
        [test_output_global, test_output_extrafields]
    )
end

"""
    copse_output_load_bergman2004(output_dir) -> test_output::PALEOmodel.OutputWriters::OutputMemory

Load data from Noam's COPSE run (header changed to remove brackets).
Convert to PALEO format.

Original unmodified data is returned in test_output["COPSE_514"].
Renormalized data is returned in test_output["global"]
"""
function copse_output_load_bergman2004(COPSE_test_output_dir)


    outputfile = joinpath(COPSE_test_output_dir, "base_514_namechange.txt")

    @info "copse_output_load_bergman2004 loading output from file $outputfile"

    COPSE_514 = DataFrames.DataFrame(CSV.File(outputfile; delim=" ", ignorerepeated=true))
    
    # Add a copy of the original data
    C_514 = Dict("COPSE_514" => COPSE_514)

    # Define normalisation values
    norm_values = (;
        P0    = 3.1000e+15,
        N0    = 4.3500e+16,
        O0    = 3.7000e+19,
        C0    = 1e+21,
        G0    = 1e+21,
        A0    = 3.1930e+18,
        PYR0  = 1e+20,
        GYP0  = 1e+20,
        S0    = 4.0000e+19,
        CAL0  = 1.3970e+19,
    )
    # define parameter values from C514
    par_values = (;
        k1_oxfrac     = 0.8600,
        k6_fepb       = 6.0000e+09,
        k7_capb       = 1.5000e+10,
        k11_landfrac  = 0.1035,
    )

    # copy into new DataFrame, renaming and renormalizing
    glb = DataFrames.DataFrame()

    glb.tmodel  = -1e6*COPSE_514.time_My_  # convert Ma to yr
    glb.O       = COPSE_514.O2*norm_values.O0
    glb.pO2PAL  = COPSE_514.O2
        
    for fld in ("A", "P", "N", "S", "C", "G", "PYR", "GYP", "CAL")
        norm0 = getfield(norm_values, Symbol(fld*"0"))
        DataFrames.insertcols!(glb, Symbol(fld) => COPSE_514[:, fld]*norm0)
    end

    glb.pCO2PAL = COPSE_514.CO2
    glb.TEMP    = COPSE_514.T .+ PB.Constants.k_CtoK
    glb.ANOX    = COPSE_514.anox
    glb.VEG     = COPSE_514.V
                
    for fld in ("silw","carbw","pyrw","gypw","mocb","locb","mccb")
        DataFrames.insertcols!(glb, Symbol(fld) => COPSE_514[:, fld]*1e12)
    end
        
    glb.mpsb    = COPSE_514.pyrb*1e12
    glb.mgsb    = COPSE_514.gypb*1e12
    glb.phosw   = COPSE_514.phsw*1e10
    glb.oxidw   = COPSE_514.oxdw*1e12
    glb.mccb_delta=COPSE_514.d13C
    glb.S_delta = COPSE_514.d34S
    # Alk ignore as not used

    # Add some derived P burial fractions
    glb.capb = par_values.k7_capb * glb.mocb/4.5e12
    glb.fepb = (par_values.k6_fepb./par_values.k1_oxfrac).*(1.0.-glb.ANOX)
    glb.mopb = glb.mocb / 250 
    glb.psea = glb.phosw .* ( 1 .- par_values.k11_landfrac.*glb.VEG)

    # Add budget checks

    glb.clc_RedoxS = 2*(glb.GYP + glb.S)
    glb.clc_RedoxC = glb.C + glb.A;
    glb.clc_RedoxNet = glb.clc_RedoxC + glb.clc_RedoxS + glb.O

    C_514["global"] = glb

    
    return PALEOmodel.OutputWriters.OutputMemory(
        [ PALEOmodel.OutputWriters.OutputMemoryDomain("global", glb) ]
    )

end

"""
    compare_copse_output(output, comparison) -> diff::DataFrame

Calculate difference between PALEO `output` and COPSE Matlab `comparison`. 

'output' is interpolated to the model times `tmodel` of `comparison`, and `diff` also uses the
`comparison` `tmodel`.


Output stats can then be calculated eg with:
    show(sort(DataFrames.describe(diff[2:end, :], :std, :min, :max, :mean), :std), allrows=true)
"""
function compare_copse_output(output, comparison)

    # output["extrafields"] = DataFrames.DataFrame()
    # output["extrafields"].tmodel = output["global"].tmodel
    # output["extrafields"].d_mccb = output["ocean"].mccb_delta - output["global"].A_delta
    # output["extrafields"].d_ocean = output["ocean"].DIC_delta - output["global"].A_delta

    compfields = [
        # outdom    field           compdom     field       op
        ("global",  "UPLIFT",       "global",   "",         :reldiff),
        ("global",  "DEGASS",       "global",   "",         :reldiff),
        ("global",  "PG",           "global",   "",         :reldiff),
        ("global",  "EVO",          "global",   "",         :reldiff),
        ("land",    "BA_AREA",      "global",   "BA",       :reldiff),
        ("global",  "GRAN",         "global",   "",         :reldiff),
        ("land",    "GRAN_AREA",    "global",   "",         :reldiff),
        ("global",  "F_EPSILON",    "global",   "",         :reldiff),
        ("global",  "COAL",         "global",   "",         :reldiff),
        ("global",  "Bforcing",     "global",   "",         :reldiff),
        ("global",  "CAL_NORM",     "global",   "",         :reldiff),
        ("global",  "ORGEVAP_AREA", "global",   "",         :reldiff),
        ("global",  "W",            "global",   "",         :reldiff),
        ("global",  "A",            "global",   "",         :reldiff),
        ("global",  "O",            "global",   "",         :reldiff),
        ("ocean",   "P",            "global",   "",         :reldiff),
        ("global",  "P",            "global",   "",         :reldiff),
        ("ocean",   "N",            "global",   "",         :reldiff),
        ("global",  "N",            "global",   "",         :reldiff),
        ("ocean",   "S",            "global",   "",         :reldiff),
        ("sedcrust","C",            "global",   "",         :reldiff),
        ("sedcrust","G",            "global",   "",         :reldiff),
        ("sedcrust","PYR",          "global",   "",         :reldiff),
        ("sedcrust","GYP",          "global",   "",         :reldiff),
        ("atm",     "pCO2PAL",      "global",   "",         :reldiff),
        ("global",  "pCO2PAL",      "global",   "",         :reldiff),
        ("atm",     "pO2PAL",       "global",   "",         :reldiff),
        ("global",  "pO2PAL",       "global",   "",         :reldiff),
        ("global",  "phi",          "global",   "",         :reldiff),
        ("global",  "TEMP",         "global",   "temp",     :reldiff),
        ("land",    "VEG",          "global",   "",         :reldiff),
        ("land",    "ignit",        "global",   "",         :diff),
        ("land",    "firef",        "global",   "",         :reldiff),
        ("land",    "oxidw",        "global",   "",         :reldiff),
        ("land",    "carbw",        "global",   "",         :reldiff),
        ("land",    "silw",         "global",   "",         :reldiff),
        ("land",    "granw",        "global",   "",         :reldiff),
        ("land",    "basw",         "global",   "",         :reldiff),
        ("land",    "phosw",        "global",   "",         :reldiff),
        ("ocean",   "newp",         "global",   "",         :reldiff),
        ("ocean",   "ANOX",         "global",   "",         :diff),
        ("global",  "ANOX",         "global",   "",         :diff),
        ("ocean",   "denit",        "global",   "",         :reldiff),
        ("ocean",   "nfix",         "global",   "",         :reldiff),
        ("ocean",   "mocb",         "global",   "",         :reldiff),
        ("fluxOceanBurial",   "flux_total_Porg",         "global",   "mopb",     :reldiff),
        ("fluxOceanBurial",   "flux_total_Pauth",         "global",   "capb",     :reldiff),
        ("fluxOceanBurial",   "flux_total_PFe",         "global",   "fepb",     :reldiff),
        ("oceanfloor","sfw_total",  "global",   "sfw",      :reldiff),        
        ("global",  "A_delta",      "global",   "delta_A",  :diff),
        ("sedcrust","C_delta",      "global",   "delta_C",  :diff),
        ("sedcrust","G_delta",      "global",   "delta_G",  :diff),
        ("ocean",   "DIC_delta",    "extrafields","delta_ocean",:diff), 
        ("ocean",   "mccb_delta",   "extrafields",   "delta_mccb",:diff),
        ("global",  "mccb_delta",   "global",   "",         :diff), 
        ("ocean",   "mocb_delta",   "extrafields","delta_mocb",:diff),
        ("ocean",   "D_oceanDIC_A", "global",   "d_ocean",  :diff),   # DIC_delta = D_oceanDIC_A + A_delta (OK)
        # ("extrafields","d_mccb",    "global",   "d_mccb",   :diff),   # mccb_delta - A_delta   (wrong)
        ("ocean",   "D_mccb_A",     "global",   "d_mccb",   :diff),   # mccb_delta - A_delta   (wrong)
        ("ocean",   "D_B_mccb_mocb","global",   "D_B",      :diff),
        ("ocean",   "S_delta",      "global",   "delta_S",  :diff),
        ("global",  "S_delta",      "global",   "",     :diff),
        ("atm",     "CO2_delta",    "extrafields","delta_atmos", :diff),
        ("land",    "locb_delta",   "extrafields","delta_locb", :diff),
        ("sedcrust", "Sr_old_ig_delta", "global",   "delta_old_ig",  :diff),
        ("ocean",   "Sr_delta",      "global",   "delta_Sr_ocean",  :diff),

    ]

    diff = DataFrames.DataFrame()

    # Find model time for comparison
    comp_tmodel = PB.get_data(comparison, "global.tmodel")
    diff.tmodel = comp_tmodel

    for (outdom, outfield, compdom, compfield, op) in compfields
        if isempty(compfield)
            compfield = outfield
        end
        if PB.has_variable(output, outdom*"."*outfield)
            if PB.has_variable(comparison, compdom*"."*compfield)
                println("comparing output $outdom.$outfield, comparison field $compdom.$compfield")
                output_tmodel = PB.get_data(output, outdom*".tmodel")
                # skip first few points if tmodel is identical, to avoid LinearInterpolation error
                firstoutputpoint = 1
                while output_tmodel[firstoutputpoint+1] == output_tmodel[firstoutputpoint]
                    firstoutputpoint += 1
                end
                firstoutputpoint == 1 || println("    starting at output index $firstoutputpoint to skip identical tmodel values")
                interp_out = Interpolations.linear_interpolation(
                    output_tmodel[firstoutputpoint:end],
                    PB.get_total.(PB.get_data(output, outdom*"."*outfield)[firstoutputpoint:end]))
                out_i = interp_out(comp_tmodel)

                comp  = PB.get_data(comparison, compdom*"."*compfield)
                if op == :reldiff
                    diffcomp = (out_i - comp)./(0.5.*(out_i+comp) .+ eps())
                elseif op == :diff
                    diffcomp = (out_i - comp)
                else
                    error("unknown op ", op)
                end
                DataFrames.insertcols!(diff, Symbol(outfield*"_"*String(op)) => diffcomp)
            else
                println("skipping comparion $outdom.$outfield, comparison field $compdom.$compfield not present")
            end
        else
            println("skipping comparion $outdom.$outfield, output field not present")
        end
    end
    
    return diff
end

"plot PALEO output and archived Matlab COPSE output"
function copse_reloaded_comparecopse(
    output, copse_output; 
    include_Sr, 
    mccb_extrafields=false,
    pager=PALEOmodel.DefaultPlotPager(),
)

    # bodge for delta_mccb which might be in two places in archived output
    mccb_var = mccb_extrafields ? "extrafields.delta_mccb" : "global.delta_mccb"

    pdef = [
        # title                         PALEO names                                             COPSE matlab names  kwargs
        ("DEGASS forcing comparison",   ["global.DEGASS"],                                      ["global.DEGASS"],   () )
        ("land area comparison",        ["land.GRAN_AREA", "land.BA_AREA"],                     "global.".*["GRAN_AREA", "BA"], ())
        ("C degass",                    ["sedcrust.ocdeg", "sedcrust.ccdeg"],                   "global.".*["ocdeg", "ccdeg"], ())
        ("S degass",                    ["sedcrust.pyrdeg", "sedcrust.gypdeg"],                 "global.".*["pyrdeg", "gypdeg"], ())
        ("S weathering",                ["land.pyrw", "land.gypw"],                             "global.".*["pyrw", "gypw"], ())
        ("S burial",                    "fluxOceanBurial.flux_total_".*["PYR", "GYP"],          "global.".*["mpsb", "mgsb"], ())
        ("P weathering",                ["land.pland", "land.psea"],                            "global.".*["pland", "psea"], ())
        ("P burial",                    "fluxOceanBurial.flux_total_".*["PFe", "Pauth", "Porg"], "global.".*["fepb", "capb", "mopb"], ())
        ("P comparison",                ["ocean.P"],                                            ["global.P"], ())
        ("N comparison",                ["ocean.N"],                                            ["global.N"], ())
        ("N cycle comparison",          ["ocean.nfix", "ocean.denit"],                          "global.".*["nfix", "denit"], ())
        ("ANOX comparison",             ["ocean.ANOX"],                                         ["global.ANOX"], ())
        ("newp comparison",             ["ocean.newp"],                                         ["global.newp"], ())
        ("Corg burial",                 ["land.locb", "ocean.mocb"],                            ["global.locb", "global.mocb"], ())
        ("pO2 comparison",              ["atm.pO2PAL"],                                         ["global.pO2PAL"], (ylabel="pO2 (PAL)",))
        ("pCO2 comparison",             ["atm.pCO2PAL"],                                        ["global.pCO2PAL"], (ylabel="pCO2 (PAL)",))
        ("Carbon isotopes comparison",  ["ocean.mccb_delta"],                                   [mccb_var], (ylabel="delta 13C (‰)",))
    ]

    pdef_Sr = [
        ("Sr old ig isotope comparison", ["sedcrust.Sr_old_ig_delta"],                          ["global.delta_old_ig"], (ylabel="d87Sr",))
        ("Sr sed comparison",           ["sedcrust.Sr_sed"],                                    ["global.Sr_sed"], (ylabel="Sr (mol)",))
        ("Sr ocean comparison",         ["ocean.Sr"],                                           ["global.Sr_ocean"], (ylabel="Sr (mol)",))
        ("Sr ocean isotope comparison", ["ocean.Sr_delta"],                                     ["global.delta_Sr_ocean"], (ylabel="d87Sr",))
    ]

    if include_Sr
        append!(pdef, pdef_Sr)
    end

    for (title, paleonames, copsenames, kwargs) in pdef
        pager(
            (Plots.plot(title=title,  output, paleonames; kwargs...);
                Plots.plot!(copse_output,    copsenames; kwargs...))
        )
    end

    pager(:newpage)  # flush output

    return nothing
end


end
