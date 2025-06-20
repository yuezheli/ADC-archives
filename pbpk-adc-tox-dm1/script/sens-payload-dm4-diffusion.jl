# date: 2/7/2025
# author: Yuezhe Li 
# purpose of this code: analysis on if decreasing payload diffusivity would result in PK-driven payload accumulation
# PK used are SAR408701, IMGN853

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)
using DataFrames

using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using Plots

include(@projectroot("model/jones_homo.jl") );
include(@projectroot("model/param-pk.jl")); 
include(@projectroot("model/init-pk.jl")); 
include(@projectroot("script/helper.jl")); 

# turn off irrelavent things
p_base.init_sR = 0.

p_sar = deepcopy(p_base);
p_sar.PS_Score = 8  # from SAR408701 fitting


p_imgn = deepcopy(p_base);
p_imgn.PS_Score = 15          # fitted from IMGN853

tspan = (-0.01, 84*DayToHour);  # hr
dose_q3w = [0, 21, 42, 63]*DayToHour

sar_q3w_dm1 = InfusionDoses(170*HT*HT/BW, dose_q3w, p_sar);
imgn_q3w_dm1 = InfusionDoses(170*HT*HT/BW, dose_q3w, p_imgn);

p_he_1 = plot(xlabel = "Time (day)", ylabel = "Payload concentration (nM)", dpi = 1000, size = (400, 400), legend = :bottom);
plot!(xticks = (0:7:21), xlims = [0, 21]);
plot!(sar_q3w_dm1.time/DayToHour, sar_q3w_dm1.pl_liver*1E3, label = "SAR408701, 170 mg/m2 Q3W", color = :turquoise3);
plot!(imgn_q3w_dm1.time/DayToHour, imgn_q3w_dm1.pl_liver*1E3, label = "IMGN853, 170 mg/m2 Q3W", color = :hotpink);
display(p_he_1)

savefig(p_he_1, @projectroot("deliv/figure/payload/sar408701-imgn853-low-payload-diffusion-he-q3w.png"));
