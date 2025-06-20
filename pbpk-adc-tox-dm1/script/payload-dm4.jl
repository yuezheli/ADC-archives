# date: 1/30/2025
# author: Yuezhe Li 
# purpose of this code: comparison of ADCs with DM4 payload 

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)
using Plots
using CSV, DataFrames, DataFramesMeta
using Statistics

include(@projectroot("model/adc-constants.jl"));

pl_tr = CSV.read(@projectroot("data/sims/homo-sar408701.csv"),DataFrame);
pl_ms = CSV.read(@projectroot("data/sims/homo-imgn853.csv"),DataFrame);

# PK between Tusamitamab ravtansine & Mirvetuximab soravtansine
pl_tr_170 = @rsubset(pl_tr, :dose .== "170 mg/m2 Q3W");
pl_tr_190 = @rsubset(pl_tr, :dose .== "190 mg/m2 Q3W");
pl_ms_170 = @rsubset(pl_ms, :dose .== "170 mg/m2 Q3W");
pl_ms_7 = @rsubset(pl_ms, :dose .== "7.0 mg/kg Q3W");

p_plasma_1 = plot(xticks = (0:7:84), xlabel = "Time (day)", ylabel = "ADC plasma concentration (nM)", dpi = 1000, size = (400, 400), legend = :bottom);
plot!(ylims = [10, 1000], yaxis = :log10);
plot!(pl_tr_170.time/DayToHour, pl_tr_170.adcplasma*1E3, label = "SAR408701, 170 mg/m2 Q3W", color = :turquoise3);
plot!(pl_ms_170.time/DayToHour, pl_ms_170.adcplasma*1E3, label = "IMGN853, 170 mg/m2 Q3W", color = :hotpink);
display(p_plasma_1)

p_he_1 = plot(xticks = (0:7:84), xlabel = "Time (day)", ylabel = "DM4 concentration (nM)", dpi = 1000, size = (400, 400), legend = :bottom);
plot!(pl_tr_170.time/DayToHour, pl_tr_170.pl_liver*1E3, label = "SAR408701, 170 mg/m2 Q3W", color = :turquoise3);
plot!(pl_ms_170.time/DayToHour, pl_ms_170.pl_liver*1E3, label = "IMGN853, 170 mg/m2 Q3W", color = :hotpink);
display(p_he_1)

savefig(p_he_1, @projectroot("deliv/figure/payload/sar408701-imgn853-he-q3w.png"));
savefig(p_plasma_1, @projectroot("deliv/figure/pk/sar408701-imgn853-q3w.png"));

p_he_2 = plot(xticks = (0:7:84), xlabel = "Time (day)", ylabel = "DM4 concentration (nM)", dpi = 1000, size = (400, 400), legend = :bottom);
plot!(pl_tr_190.time/DayToHour, pl_tr_190.pl_liver*1E3, label = "SAR408701, 190 mg/m2 Q3W", color = :turquoise3);
plot!(pl_ms_7.time/DayToHour, pl_ms_7.pl_liver*1E3, label = "IMGN853, 7 mg/kg Q3W", color = :blue);
display(p_he_2)

savefig(p_he_2, @projectroot("deliv/figure/payload/sar408701-imgn853-he-q3w-2.png"));

# payload between 2 MTD of Tusamitamab ravtansine
df2 = vcat(
    @rsubset(pl_tr, :dose .== "100 mg/m2 Q2W"), 
    @rsubset(pl_tr, :dose .== "170 mg/m2 Q3W"), 
);

plot2_pl_he = plot(xlabel = "Time (day)", ylabel = "DM4 concentration (nM)", dpi = 1000, size = (400, 400), ylims = [0, 3.5], xticks = 0:7:84 ); 
plot!(legend = :top, legendcolumns=2, legendfontsize = 8, xlabel = "Time (day)", ylabel = "DM4 conc (nM)");
plot!(df2.time/DayToHour, df2.pl_liver*1E3, group = df2.dose, palette = palette([:purple, :green], 2), alpha = 0.8);
display(plot2_pl_he)

savefig(plot2_pl_he, @projectroot("deliv/figure/payload/sar408701-he-q2w-q3w.png"));

# summary statistics comparison 

Cavg = [mean(@rsubset(pl_tr, :dose .== "100 mg/m2 Q2W").pl_liver), 
        mean(@rsubset(pl_tr, :dose .== "170 mg/m2 Q3W").pl_liver), 
        mean(@rsubset(pl_ms, :dose .== "170 mg/m2 Q3W").pl_liver), 
        mean(@rsubset(pl_ms, :dose .== "6.0 mg/kg Q3W").pl_liver), 
        mean(@rsubset(pl_ms, :dose .== "7.0 mg/kg Q3W").pl_liver), 
        mean(@rsubset(pl_tr, :dose .== "190 mg/m2 Q3W").pl_liver), ]

Cmax = [maximum(@rsubset(pl_tr, :dose .== "100 mg/m2 Q2W").pl_liver), 
        maximum(@rsubset(pl_tr, :dose .== "170 mg/m2 Q3W").pl_liver), 
        maximum(@rsubset(pl_ms, :dose .== "170 mg/m2 Q3W").pl_liver), 
        maximum(@rsubset(pl_ms, :dose .== "6.0 mg/kg Q3W").pl_liver), 
        maximum(@rsubset(pl_ms, :dose .== "7.0 mg/kg Q3W").pl_liver), 
        maximum(@rsubset(pl_tr, :dose .== "190 mg/m2 Q3W").pl_liver)]
