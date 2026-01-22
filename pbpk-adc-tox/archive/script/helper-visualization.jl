# date: 6/30/2025 
# author: Yuezhe Li 
# purpose of this code: to visualize data

# plot plasma concentration in 
function PlotSimulationPlasma(dict_sim, dict_obs; adc_name = missing, colorPALETTE = :batlow10, xrange = missing, yrange = missing, ylog = false, legendcolumnsnum = 1)
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
        plot!(dict_sim[dose_i].t/hr_per_day, [dict_sim[dose_i].u[i].C_EXG_Plasma for i in 1:length(dict_sim[dose_i].t)], label="sims, $(round(dose_i, digits=2)) mgkg", width=2, color = usepalette[i], aplha = 0.7)
        plot!(dict_obs[dose_i].time_d, dict_obs[dose_i].ADC_uM, label="obs, $(round(dose_i, digits=2)) mgkg", seriestype = :scatter, color = usepalette[i], alpha = 0.6); 
    end
    return plt__
end


# plot payload concentration in endothelial cells
function PlotSimulationLiverEndoCyto(dict_sim; adc_name = missing, colorPALETTE = :batlow10, xrange = missing, yrange = missing, ylog = false, pl_ic50 = missing)
    ylabel_text = "Liver endothelial " * adc_name * " (uM)"
    dose_mgkg_list = collect(keys(dict_sim))
    usepalette = cgrad(colorPALETTE, length(dose_mgkg_list))
    plt__ = plot(xlabel = "Time (d)", ylabel = ylabel_text, dpi = 300, size = (400,400), background_color_legend = nothing);
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
        plot!(dict_sim[dose_i].t/hr_per_day, [dict_sim[dose_i].u[i].end_cyto_payload[2] for i in 1:length(dict_sim[dose_i].t)], label="sims, $(round(dose_i, digits=2)) mgkg", width=2, color = usepalette[i], aplha = 0.7)
    end
    if !ismissing(pl_ic50)
        hline!([pl_ic50], label = "PL IC50", linestyle = :dashdot, color = "black", alpha = 0.4); 
    end
    return plt__
end
