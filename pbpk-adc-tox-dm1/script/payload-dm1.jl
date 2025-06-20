# date: 1/28/2025
# author: Yuezhe Li 

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)
using Plots
using CSV, DataFrames, DataFramesMeta
using Statistics

include(@projectroot("model/adc-constants.jl"));

pl_cm = CSV.read(@projectroot("data/sims/homo-cantuzumab-mertansine.csv"),DataFrame);
pl_tdm1 = CSV.read(@projectroot("data/sims/homo-tdm1.csv"),DataFrame);
pl_pca062 = CSV.read(@projectroot("data/sims/homo-pca062.csv"),DataFrame);
pl_cm = CSV.read(@projectroot("data/sims/homo-cantuzumab-mertansine.csv"),DataFrame);

# CM only 
Cmax = [maximum(@rsubset(pl_cm, :dose == "235 mg/m2 Q3W").pl_liver), 
        maximum(@rsubset(pl_cm, :dose == "115 mg/m2 QW").pl_liver), 
        maximum(@rsubset(pl_cm, :dose == "45 mg/m2 QOD").pl_liver)]

Cavg = [mean(@rsubset(pl_cm, :dose == "235 mg/m2 Q3W").pl_liver), 
        mean(@rsubset(pl_cm, :dose == "115 mg/m2 QW").pl_liver), 
        mean(@rsubset(pl_cm, :dose == "45 mg/m2 QOD").pl_liver)]


# payload conc 
cm235 = @rsubset(pl_cm, :dose == "235 mg/m2 Q3W");
cm115 = @rsubset(pl_cm, :dose == "115 mg/m2 QW");
cm45 = @rsubset(pl_cm, :dose == "45 mg/m2 QOD");

pliver1 = plot(xticks = (0:7:84), xlabel = "Time (day)", ylabel = "DM1 concentration (nM)", dpi = 1000, size = (400, 400), legend = :topleft, ylims = [0, 50])
plot!(cm235.time/DayToHour, cm235.pl_liver*1E3, label = "235 mg/m2 Q3W");
plot!(cm115.time/DayToHour, cm115.pl_liver*1E3, label = "115 mg/m2 QW");
plot!(cm45.time/DayToHour, cm45.pl_liver*1E3, label = "45 mg/m2 3 times a week, 3 in 4 weeks");
pliver1

savefig(pliver1, @projectroot("deliv/figure/payload/cm-235q3w-115qw-liver.png"));

# across DM1 comparison 
cm3point6 = @rsubset(pl_cm, :dose == "3.6 mg/kg Q3W");
@rsubset!(cm3point6, :time <= 42*DayToHour);
tdm1_3point6 = @rsubset(pl_tdm1, :dose == "3.6mg/kg");
@rsubset!(tdm1_3point6, :time <= 42*DayToHour);

Cmax = [maximum(@rsubset(pl_tdm1, :dose == "3.6mg/kg").pl_liver), 
        maximum(@rsubset(pl_cm, :dose == "3.6 mg/kg Q3W").pl_liver)]


Cavg = [mean(@rsubset(pl_tdm1, :dose == "3.6mg/kg").pl_liver), 
        mean(@rsubset(pl_cm, :dose == "3.6 mg/kg Q3W").pl_liver)]

pliver2 = plot(xticks = (0:7:42), xlabel = "Time (day)", ylabel = "DM1 concentration (nM)", dpi = 1000, size = (400, 400), legend = :topleft, ylims = [0, 15]);
plot!(palette = palette([:turquoise3, :hotpink], 2));
plot!(cm3point6.time/DayToHour, cm3point6.pl_liver*1E3, label = "CM, 3.6 mg/kg Q3W");
plot!(tdm1_3point6.time/DayToHour, tdm1_3point6.pl_liver*1E3, label = "T-DM1, 3.6 mg/kg Q3W");
display(pliver2)

savefig(pliver2, @projectroot("deliv/figure/payload/cm-tdm1-liver.png"));

pplasma2 = plot(xticks = (0:7:42), xlabel = "Time (day)", ylabel = "ADC plasma concentration (nM)", dpi = 1000, size = (400, 400), legend = :topright, ylims = [0, 450]);
plot!(palette = palette([:turquoise3, :hotpink], 2));
plot!(cm3point6.time/DayToHour, cm3point6.adcplasma*1E3, label = "CM, 3.6 mg/kg Q3W");
plot!(tdm1_3point6.time/DayToHour, tdm1_3point6.adcplasma*1E3, label = "T-DM1, 3.6 mg/kg Q3W");
display(pplasma2)

savefig(pplasma2, @projectroot("deliv/figure/pk/cm-tdm1.png"));

