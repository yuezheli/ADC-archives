# date: 1/29/2025
# author: Yuezhe Li 
# purpose of this code: to plot payload concentration of T-DM1 and T-Dxd

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)
using Plots
using CSV, DataFrames, DataFramesMeta
using Statistics

include(@projectroot("model/adc-constants.jl"));

pl_tdxd = CSV.read(@projectroot("data/sims/homo-tdxd.csv"),DataFrame);
pl_tdm1 = CSV.read(@projectroot("data/sims/homo-tdm1.csv"),DataFrame);

@rsubset!(pl_tdxd, :time .<= 42*DayToHour);

pl_tdxd_5point4 = @rsubset(pl_tdxd, :dose == "5.4mg/kg")
pl_tdm1_3point6 = @rsubset(pl_tdm1, :dose == "3.6mg/kg")

# liver tox 
pliver1 = plot(xlabel = "Time (day)", ylabel = "Payload concentration (nM)", dpi = 1000, size = (400, 400), legend = :bottomright);
plot!(xticks = (0:7:42), ylims = [0.005, 10], yaxis = :log, yticks = [0.01, 0.1, 1, 10]);
plot!(pl_tdxd_5point4.time/DayToHour, pl_tdxd_5point4.pl_liver*1E3, label = "T-Dxd 5.4 mg/kg Q3W", color = :deepskyblue2);
plot!(pl_tdm1_3point6.time/DayToHour, pl_tdm1_3point6.pl_liver*1E3, label = "T-DM1 3.6 mg/kg Q3W", color = :mediumpurple);
hline!([3], linestyle=:dashdot, label="IC50, DM1", color = :mediumpurple, alpha = 0.7); # https://www.mdpi.com/2218-273X/10/3/360
hline!([0.371], linestyle=:dashdot, label="IC50, Dxd", color = :deepskyblue2, alpha = 0.7); # https://aacrjournals.org/mct/article/21/4/635/689570/DS-7300a-a-DNA-Topoisomerase-I-Inhibitor-DXd-Based
display(pliver1)

savefig(pliver1, @projectroot("deliv/figure/payload/tmd1-tdxd-liver.png"));

# compute Cmax and Cavg 
Cmax = [maximum(pl_tdxd_5point4.pl_liver)*1E3, maximum(pl_tdm1_3point6.pl_liver)*1E3]

Cavg = [mean(pl_tdxd_5point4.pl_liver)*1E3, mean(pl_tdm1_3point6.pl_liver)*1E3]

# skin tox 
pskin1 = plot(xlabel = "Time (day)", ylabel = "Payload concentration (nM)", dpi = 1000, size = (400, 400), legend = :topright);
plot!(xticks = (0:7:42), ylims = [0.005, 10], yaxis = :log, yticks = [0.01, 0.1, 1, 10]);
plot!(pl_tdxd_5point4.time/DayToHour, pl_tdxd_5point4.pl_skin*1E3, label = "T-Dxd 5.4 mg/kg Q3W", color = :deepskyblue2);
plot!(pl_tdm1_3point6.time/DayToHour, pl_tdm1_3point6.pl_skin*1E3, label = "T-DM1 3.6 mg/kg Q3W", color = :mediumpurple);
hline!([3], linestyle=:dashdot, label="IC50, DM1", color = :mediumpurple, alpha = 0.7); # https://www.mdpi.com/2218-273X/10/3/360
hline!([0.371], linestyle=:dashdot, label="IC50, Dxd", color = :deepskyblue2, alpha = 0.7); # https://aacrjournals.org/mct/article/21/4/635/689570/DS-7300a-a-DNA-Topoisomerase-I-Inhibitor-DXd-Based
display(pskin1)

savefig(pskin1, @projectroot("deliv/figure/payload/tmd1-tdxd-skin.png"));

Cmax = [maximum(pl_tdxd_5point4.pl_skin)*1E3, maximum(pl_tdm1_3point6.pl_skin)*1E3]

Cavg = [mean(pl_tdxd_5point4.pl_skin)*1E3, mean(pl_tdm1_3point6.pl_skin)*1E3]
