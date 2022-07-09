
import PALEOboxes as PB

import PALEOmodel


function copse_reloaded_reloaded_expts(
    basemodel, expts;
    modelpars=Dict{}(),
)

    ####################
    # set basemodel
    ####################
    if basemodel == "reloaded"
        
        # modelpars = Dict("CIsotope"=>"ScalarData", "SIsotope"=>"ScalarData", "CIsotopeReacts"=>false)       
        model = PB.create_model_from_config(
            joinpath(@__DIR__, "COPSE_reloaded_reloaded_cfg.yaml"), 
            "model1",
            modelpars=modelpars,
            # sort_methods_algorithm=PB.dfs_methods
        )
        
    else
        error("unknown basemodel ", basemodel)
    end


    ###############################################
    # choose an 'expt' (a delta to the base model)
    ###############################################
    for expt in expts
        if expt == "baseline"
            # baseline configuration

        elseif length(expt)==2 && expt[1] == "tforce_constant"
            # Set constant forcing time
            # NB: use with 'modelpars'=Dict("tforcevar"=>"tforce_constant")        
            PB.set_variable_attribute!(model, "global", "tforce_constant", :initial_value, expt[2])
            
        # Experiments in Figs. 10-12 of COPSE reloaded (for basemodel 'reloaded')   
        elseif  expt == "climsens15"
            # Fig.11 blue
            PB.set_parameter_value!(model, "global", "temp_global", "k_c", 2.164) # 1.5        
        elseif  expt == "climsens225"
            # Fig.11 green
            PB.set_parameter_value!(model, "global", "temp_global", "k_c", 3.246) # 2.25        
        elseif  expt == "climsens45"
            # Fig.11 cyan
            PB.set_parameter_value!(model, "global", "temp_global", "k_c", 6.492) # 4.5        
        elseif  expt == "climsens6"
            # Fig.11 red
            PB.set_parameter_value!(model, "global", "temp_global", "k_c", 8.656) # 6

        elseif  expt == "VCI"
            # Fig.12 magenta
            PB.set_parameter_value!(model, "ocean", "oceanburial_copse", "f_CPsea", "VCI")

        elseif  expt == "mocbProdLinear"
            # not shown in paper - mocb linearly proportional to new production
            PB.set_parameter_value!(model, "ocean", "oceanburial_copse", "f_mocb_b", 1.0)

        elseif  expt == "noNcycle"
            # not shown in paper
            PB.set_parameter_value!(model, "ocean", "marinebiota_copse", "f_ncycle", false)
            PB.set_parameter_value!(model, "ocean", "oceanburial_copse", "f_ncycle", false)

        elseif length(expt) == 5 && expt[1] == "setpar"
            # generic parameter set (setpar, <domain>, <reaction>, <parname>, <parvalue)
            _, domname, reactname, parname, parvalue = expt            
            PB.set_parameter_value!(model, domname, reactname, parname, parvalue)

        elseif length(expt)==4 && expt[1] == "CO2pulse"
            _, size, pstart, duration = expt # (_, mol C, yr, yr)
            
            CO2pulse = PB.get_reaction(model, "global", "CO2pulse")
            delta = -5       
            
            @info "expt $expt apply CO2pulse size $size (mol C) delta $delta (per mil) duration $duration (yr) start time $pstart (model yr)"
            # Witches hat perturbation
            PB.set_parameter_value!(CO2pulse, "perturb_times",  
                [-1e30, pstart, pstart+duration/2, pstart+duration,  1e30]) 
            PB.set_parameter_value!(CO2pulse, "perturb_totals", 
                [0.0,   0.0,    2*size/duration,    0.0,             0.0])
            PB.set_parameter_value!(CO2pulse, "perturb_deltas",
                delta.*[1.0,   1.0,    1.0,         1.0,             1.0])
    

        else
            error("unknown expt ", expt)
        end
    end

    # set parameters for carbonate-silicate weathering steady-state
    rct_ocean       = PB.get_reaction(model, "ocean", "oceanburial_copse")
    mtotpb = rct_ocean.pars.k2_mocb.v/rct_ocean.pars.CPsea0.v  +rct_ocean.pars.k7_capb.v +rct_ocean.pars.k6_fepb.v
    rct_degass      = PB.get_reaction(model, "sedcrust", "sedcrust_copse")
    rct_sfw         = PB.get_reaction(model, "oceanfloor", "sfw")
    rct_land_fluxes = PB.get_reaction(model, "land", "land_fluxes")
    PALEOcopse.Land.LandCOPSEReloaded.set_steady_state(
        rct_land_fluxes, 
        rct_ocean.pars.k2_mocb.v,
        mtotpb,
        rct_degass.pars.k12_ccdeg.v,
        rct_degass.pars.k13_ocdeg.v, 
        0.0,
        rct_sfw.pars.k_sfw.v
    )

    # set parameters for Sr steady state
    rct_Sr_land       = PB.get_reaction(model, "land", "Sr_land")
    rct_Sr_oceanfloor = PB.get_reaction(model, "oceanfloor", "Sr_oceanfloor")
    PALEOcopse.BioGeoChem.Strontium.set_Sr_fluxes_steady_state!(
        rct_Sr_land,
        rct_Sr_oceanfloor,
        rct_land_fluxes.pars.k_basw.v,
        rct_land_fluxes.pars.k_granw.v, 
        rct_land_fluxes.pars.k14_carbw.v, 
        rct_sfw.pars.k_sfw.v
    )
                                          
    run = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

    return run
end

function copse_reloaded_reloaded_plot_summary(
    outputs; 
    pager=PALEOmodel.DefaultPlotPager(),
    extrakwargs::NamedTuple=NamedTuple()
)

    pager(
        plot(title="O2", 100*PALEOmodel.get_array.(outputs, "land.mrO2"); ylabel="O_2 (%)", extrakwargs...),

        plot(title="pCO2", 1e6*PALEOmodel.get_array.(outputs, "atm.pCO2atm"); ylabel="pCO2 (ppm)", extrakwargs...),

        # ocean volume in litres = PB.Constants.k18_oceanmass/1.027 = 1.397e21 / 1.027 = 1.360e21
        # [SO4] mM = 1e3 * total SO4 (mol) / ocean vol (litres) = 7.351e-19 * SO4 (mol)
        plot(title="[SO4]", 7.351e-19*PALEOmodel.get_array.(outputs, "ocean.S"); ylabel="[SO4] mM", extrakwargs...),

        plot(title="Carbon isotopes",  outputs, "ocean.mccb_delta"; ylabel="delta 13C (‰)", extrakwargs...),
           
        plot(title="d34S SO4", outputs, "ocean.S_delta"; ylabel="d34S (‰)", extrakwargs...),
  
        plot(title="d87/d86 Sr", outputs, ["ocean.Sr_delta"]; ylabel="d87/d86 Sr", extrakwargs...),
    
        :newpage,
    )
end

function copse_reloaded_reloaded_plot(
    output; 
    pager=PALEOmodel.DefaultPlotPager(),
    extrakwargs::NamedTuple=NamedTuple()
)
    if isa(output, AbstractVector)
        have_CO2pulse = all(PB.has_variable.(output, "global.CO2pulse"))
    else
        have_CO2pulse = PB.has_variable(output, "global.CO2pulse")
    end

    # Conservation checks
    pager(
        plot(title="Total carbon",          output, ["global.total_C"], ylabel="mol C"; extrakwargs...),
        plot(title="Total carbon moldelta", output, ["global.total_C.v_moldelta"], ylabel="mol C * per mil"; extrakwargs...),
        plot(title="Total sulphur",         output, ["global.total_S"], ylabel="mol S"; extrakwargs...),
        plot(title="Total sulphur moldelta",output, ["global.total_S.v_moldelta"], ylabel="mol S * per mil"; extrakwargs...),
        plot(title="Total redox",           output, ["global.total_redox"], ylabel="mol O2 equiv"; extrakwargs...),
        :newpage,

        # Forcings
        plot(title="Physical forcings",     output, "global.".*["DEGASS", "UPLIFT", "PG"], ylim=(0, 2.0), ylabel="normalized forcing"; extrakwargs...),
        plot(title="Land area forcings",    output, ["land.BA_AREA", "land.GRAN_AREA", "global.ORGEVAP_AREA"],  ylim=(0, 2.5), ylabel="normalized forcing"; extrakwargs...),
        plot(title="Basalt area forcings",  output, ["global.CFB_area", "land.oib_area"],  ylabel="area (km^2)",),
        plot(title="Evolutionary forcings", output, "global.".*["EVO", "W", "Bforcing", "CPland_relative", "F_EPSILON", "COAL"], ylabel="normalized forcing"; extrakwargs...),    
    
        have_CO2pulse ?
            plot(title="CO2pulse", output,  ["global.CO2pulse"],  ylabel="CO2 pulse (mol C yr-1)"; extrakwargs...) : :skip,

        # Outputs
        plot(title="pCO2",                  output, "atm.".*["pCO2PAL"],  ylabel="pCO2 (PAL)"; extrakwargs...),
        plot(title="Temperature",           output, ["global.TEMP"],  ylabel="T (K)"; extrakwargs...),
        plot(title="Oxygen",                output, ["atm.pO2PAL", "ocean.ANOX"]; extrakwargs...),
        plot(title="Carbon",                output, ["global.total_C", "sedcrust.C", "sedcrust.G", "atmocean.A"], ylabel="mol C"; extrakwargs...),
        plot(title="Sulphur",               output, ["global.total_S", "ocean.S", "sedcrust.PYR", "sedcrust.GYP"], ylabel="mol S"; extrakwargs...),
        plot(title="Silicate weathering",   output, ["land.silw", "land.granw", "land.basw", "oceanfloor.sfw_total"], ylabel="flux (mol/yr)"; extrakwargs...),
        plot(title="Phosphorus weathering", output, ["land.phosw", "land.phosw_c", "land.phosw_o", "land.phosw_s"], ylabel="flux (mol P/yr)"; extrakwargs...),
        plot(title="P ocean reservoir",     output, ["ocean.P"], ylabel="P (mol)"; extrakwargs...), 
        plot(title="organic C burial",      output, ["ocean.mocb", "land.locb"], ylabel="Corg burial (mol C yr-1)"; extrakwargs...), 
        plot(title="S burial",              output, ["fluxOceanBurial.flux_total_GYP", "fluxOceanBurial.flux_total_PYR"], ylabel="S burial (mol S yr-1)"; extrakwargs...), 
        plot(title="Land biota",            output, ["land.V_T", "land.V_o2", "land.V_co2", "land.V_npp", "land.firef", "land.VEG", "global.COAL", "atm.pO2PAL"]; extrakwargs...),
        plot(title="Sulphur isotopes",      output, ["ocean.S_delta", "sedcrust.PYR_delta", "sedcrust.GYP_delta"], ylabel="delta 34S (per mil)"; extrakwargs...),
        plot(title="Carbon isotopes",       output, ["atmocean.A_delta", "ocean.DIC_delta", "atm.CO2_delta", "ocean.mccb_delta", "sedcrust.C_delta"], ylabel="delta 13C (per mil)"; extrakwargs...),
        plot(title="Sr sed reservoir",      output, ["sedcrust.Sr_sed"], ylabel="Sr_sed (mol)"; extrakwargs...),    
        plot(title="Sr ocean reservoir",    output, ["ocean.Sr"], ylabel="Sr (mol)"; extrakwargs...), 
        plot(title="Sr fluxes",             output, ["fluxRtoOcean.flux_Sr", "fluxOceanBurial.flux_total_Sr", "fluxOceanfloor.soluteflux_Sr" ], ylabel="Sr flux (mol yr-1)"; extrakwargs...),
        plot(title="Sr isotopes",           output, ["sedcrust.Sr_mantle_delta", "sedcrust.Sr_new_ig_delta", "sedcrust.Sr_old_ig_delta", "sedcrust.Sr_sed_delta", "ocean.Sr_delta"], ylabel="87Sr"; extrakwargs...),

        :newpage,  # flush partial screen
    )
    return nothing
end


