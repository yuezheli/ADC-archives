# date: 1/27/2025
# author: Yuezhe Li 
# purpose of this code: analysis for PBD-based payload 

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)
using Plots
using CSV, DataFrames, DataFramesMeta
using Statistics

include(@projectroot("script/helper.jl")); 

pl_medi2228 = CSV.read(@projectroot("data/sims/homo-medi2228.csv"),DataFrame);
pl_adct601 = CSV.read(@projectroot("data/sims/homo-adct601.csv"),DataFrame);
pl_adct402 = CSV.read(@projectroot("data/sims/homo-adct402.csv"),DataFrame);

# 200 ug/kg dose (15 mg for ADCT-601)
pl_2228_200 = @rsubset(pl_medi2228, :dose == "0.2 mg/kg");
pl_601_200 = @rsubset(pl_adct601, :dose == "15.0 mg");
pl_402_200 = @rsubset(pl_adct402, :dose == "200 ug/kg");

pskin_200 = plot(xticks = 0:7:42, title = "PBD ADC dosed at 0.2 mg/kg Q3W", titlefontsize = 10, legend = :topleft, legendcolumns=1, dpi = 1000, size = (400, 400)); 
plot!(xlabel = "Time (day)", ylabel = "skin endothelial payload concentration (pM)", ylims = [0, 1], guidefontsize = 10);
plot!(pl_2228_200.time/DayToHour, pl_2228_200.pl_skin*1E6, color = RGB(41/255, 127/255, 157/255), label = "MEDI2228");
plot!(pl_601_200.time/DayToHour, pl_601_200.pl_skin*1E6, color = RGB(58/255, 56/255, 56/255), label = "Mipasetamab Uzoptirine");
plot!(pl_402_200.time/DayToHour, pl_402_200.pl_skin*1E6, color = RGB(255/255, 123/255, 118/255), label = "Loncastuximab Tesirine");
display(pskin_200)

pliver_200 = plot(xticks = 0:7:42, title = "PBD ADC dosed at 0.2 mg/kg Q3W", titlefontsize = 10, legend = :topleft, legendcolumns=1, dpi = 1000, size = (400, 400)); 
plot!(xlabel = "Time (day)", ylabel = "liver endothelial payload concentration (pM)", ylims = [0, 5], guidefontsize = 10);
plot!(pl_2228_200.time/DayToHour, pl_2228_200.pl_liver*1E6, color = RGB(41/255, 127/255, 157/255), label = "MEDI2228");
plot!(pl_601_200.time/DayToHour, pl_601_200.pl_liver*1E6, color = RGB(58/255, 56/255, 56/255), label = "Mipasetamab Uzoptirine");
plot!(pl_402_200.time/DayToHour, pl_402_200.pl_liver*1E6, color = RGB(255/255, 123/255, 118/255), label = "Loncastuximab Tesirine");
display(pliver_200)

ctrough_he = [Ctrough(pl_2228_200, [42], "liver"), Ctrough(pl_601_200, [42], "liver"), Ctrough(pl_402_200, [42], "liver")] 
cmax_he = [maximum(pl_2228_200.pl_liver), maximum(pl_601_200.pl_liver), maximum(pl_402_200.pl_liver)]
cavg_he = [mean(pl_2228_200.pl_liver), mean(pl_601_200.pl_liver), mean(pl_402_200.pl_liver)]

savefig(pskin_200, @projectroot("deliv/figure/payload/pbd-payload-q3w-skin-200.png"));
savefig(pliver_200, @projectroot("deliv/figure/payload/pbd-payload-q3w-liver-200.png"));
