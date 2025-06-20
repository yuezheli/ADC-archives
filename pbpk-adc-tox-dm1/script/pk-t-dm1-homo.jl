# date: 1/6/2025
# author: Yuezhe Li 
# purpose of this code: to simulate for T-DM1 human PK 

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)
using DataFrames, CSV
# read observed 
pk_yamamoto =  CSV.read("data/yamamoto-2015.csv",DataFrame);
pk_girish =  CSV.read("data/girish-2012.csv",DataFrame);

pk_tdm1 = vcat(pk_yamamoto, pk_girish)

using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using Plots

# simulation 
include(@projectroot("model/jones_homo.jl") );
include(@projectroot("model/param-pk.jl")); 
include(@projectroot("model/init-pk.jl")); 
include(@projectroot("script/helper.jl"));

p_tdm1 = deepcopy(p_base);
p_tdm1.PS_Score = 6.
p_tdm1.init_sR = 0.004

tspan = (-0.01, DayToHour*42);      # [hr]
AddDose = [0., 21.] * DayToHour  # [hr]

first_dose = [0.3, 0.6, 1.2, 2.4, 3.6, 4.8, 1.8]; # [mg/kg]

df_tdm1 = DataFrame(time = [], pl_lung = [], pl_liver = [], pl_heart = [], pl_muscle = [], pl_skin = [], pl_adipose = [], pl_bone = [], pl_brain = [], 
                  pl_kidney = [], pl_si = [], pl_li = [], pl_pancreas = [], pl_thymus = [], pl_spleen = [], adcplasma = [], dose = []);

for init_dose in first_dose
    tmpdf = InfusionDoses(init_dose, AddDose, p_tdm1, infusion_time = 0.25);
    tmpdf.dose .= string(init_dose) * "mg/kg"
    df_tdm1 = vcat(df_tdm1, tmpdf)
end

# PK

p_pk_tdm1 = plot(legend = :topright, ylims = [1E-1, 1E3], yaxis = :log, size=(350,350), legendcolumns=3, palette = :Paired_7, dpi = 1000);
plot!(df_tdm1.time/DayToHour, df_tdm1.adcplasma*MW/1E3, group = df_tdm1.dose, linewidth = 2, alpha = 0.8); 
scatter!(pk_tdm1.time_day, pk_tdm1.T_DM1_ugperml, group = pk_tdm1.Dose, ma = 0.6, label = false, linewidth=0);
xlabel!("Time (day)"); xlims!(0, 30); xticks!(0:7:63); ylabel!("plasma T-DM1 concentration (ug/mL)"); 
display(p_pk_tdm1)

savefig(p_pk_tdm1, @projectroot("deliv/figure/pk/t-dm1-homo.png"));

# PAYLOAD 

CSV.write(@projectroot("data/sims/homo-tdm1.csv"), df_tdm1)

plung = plot(df_tdm1.time/DayToHour, df_tdm1.pl_lung, group = df_tdm1.dose, title = "lung endothelial payload, tmd1")
pliver = plot(df_tdm1.time/DayToHour, df_tdm1.pl_liver, group = df_tdm1.dose, title = "liver endothelial payload, T-DM1", xticks = (0:7:63))

savefig(pliver, @projectroot("deliv/figure/payload/t-dm1-q3w-liver.png"));
