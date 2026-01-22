# date: 1/13/2026 
# author: Yuezhe Li 
# purpose of this code: script to host visualization functions

function PlotSimulationPlasma(dict_sim, dict_obs, mdl;
    adc_name = missing, colorPALETTE = :batlow10, xrange = missing, yrange = missing, ylog = false, legendcolumnsnum = 1)
    
    ylabel_text = "Plasma " * adc_name * " (uM)"
    dose_mgkg_list = collect(keys(dict_sim))
    usepalette = cgrad(colorPALETTE, length(dose_mgkg_list))

    plt__ = plot(xlabel = "Time (d)", ylabel = ylabel_text, dpi = 300, size = (400,400), background_color_legend = nothing, legend_columns = legendcolumnsnum);
    if !ismissing(xrange)
        plot!(xlims = xrange);
    end
    if !ismissing(yrange)
        plot!(ylims = yrange);
    end
    if ylog
        plot!(yaxis = :log);
    end
    for (i, dose_i) in enumerate(dose_mgkg_list)
        plot!(dict_sim[dose_i].t/hr_per_day, dict_sim[dose_i][mdl.plasma_exg.C_Plasma], label="sims, $(round(dose_i, digits=2)) mgkg", width=2, color = usepalette[i], aplha = 0.7)
        plot!(dict_obs[dose_i].time_d, dict_obs[dose_i].ADC_uM, label="obs, $(round(dose_i, digits=2)) mgkg", seriestype = :scatter, color = usepalette[i], alpha = 0.6); 
    end
    plot!(xticks = [0, 7, 14, 21, 28, 35, 42, 49, 56, 63]);
    return plt__
end

function PlotSimulationPlasmaPL(dict_sim, dict_obs, mdl;
    adc_name = missing, colorPALETTE = :batlow10, xrange = missing, yrange = missing, ylog = false, legendcolumnsnum = 1)
    
    ylabel_text = "Plasma " * adc_name * " (uM)"
    dose_mgkg_list = collect(keys(dict_sim))
    usepalette = cgrad(colorPALETTE, length(dose_mgkg_list))
    plt__ = plot(xlabel = "Time (d)", ylabel = ylabel_text, dpi = 300, size = (400,400), background_color_legend = nothing, legend_columns = legendcolumnsnum);
    if !ismissing(xrange)
        plot!(xlims = xrange);
    end
    if !ismissing(yrange)
        plot!(ylims = yrange);
    end
    if ylog
        plot!(yaxis = :log);
    end
    for (i, dose_i) in enumerate(dose_mgkg_list)
        plot!(dict_sim[dose_i].t/hr_per_day, dict_sim[dose_i][mdl.plasma_pl.C_PL_Plasma], label="sims, $(round(dose_i, digits=2)) mgkg", width=2, color = usepalette[i], aplha = 0.7)
        plot!(dict_obs[dose_i].time_d, dict_obs[dose_i].pl_uM, label="obs, $(round(dose_i, digits=2)) mgkg", seriestype = :scatter, color = usepalette[i], alpha = 0.6); 
    end
    plot!(xticks = [0, 7, 14, 21, 28, 35, 42, 49, 56, 63]);
    return plt__
end
