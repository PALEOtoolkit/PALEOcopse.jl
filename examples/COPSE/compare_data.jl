include("../../DataCompilations/DataCompilations.jl")

"plot PALEO output and comparison datasets"
function copse_reloaded_comparedata(
    outputs::Vector{<: PALEOmodel.AbstractOutputWriter}; 
    include_Sr, 
    pager=PALEOmodel.DefaultPlotPager(),
)
    
    pager(
        (plot(title="O2 charcoal", 100*PALEOmodel.get_array.(outputs, "land.mrO2"), ylabel="O_2 (%)",);
            DataCompilations.plot_O2_charcoal()),

        (plot(title="pCO2 Royer etal (2014)", 1e6*PALEOmodel.get_array.(outputs, "atm.pCO2atm"), ylabel="pCO2 (ppm)",);
            DataCompilations.plot_royer_CO2()),

        # ocean volume in litres = PB.Constants.k18_oceanmass/1.027 = 1.397e21 / 1.027 = 1.360e21
        # [SO4] mM = 1e3 * total SO4 (mol) / ocean vol (litres) = 7.351e-19 * SO4 (mol)
        (plot(title="[SO4] fluid inclusion Algeo etal (2015)", 7.351e-19*PALEOmodel.get_array.(outputs, "ocean.S"), ylabel="[SO4] mM",);
            DataCompilations.plot_SO4()),

        (plot(title="Carbon isotopes BM2017",  outputs, "ocean.mccb_delta", ylabel="delta 13C (‰)",);
            DataCompilations.plot_d13C_BM2017()),

        (plot(title="Carbon isotopes Saltzman", outputs, "ocean.mccb_delta", ylabel="delta 13C (‰)",);
            DataCompilations.plot_d13C_Saltzman()),

        (plot(title="d34S Algeo etal (2015)", outputs, "ocean.S_delta", ylabel="d34S (‰)",);
            DataCompilations.plot_d34s_algeo2015()),

        (plot(title="d34S Algeo etal (2015) 9-pt mov av", outputs, ["ocean.S_delta"], ylabel="d34S (‰)",);
            DataCompilations.plot_d34s_algeo2015_mov(k=9)),
    )

    if include_Sr
        pager(
            (plot(title="d87/d86 Sr LOWESS, Cox etal (2016)", outputs, ["ocean.Sr_delta"], ylabel="d87/d86 Sr",);
                DataCompilations.plot_d8786Sr_LOWESS();
                DataCompilations.plot_d8786Sr_cox_mov(k=9)),
        )
    end
    
    pager(:newpage) # flush output

    return nothing
end


