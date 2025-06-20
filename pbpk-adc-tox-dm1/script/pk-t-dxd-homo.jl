# date: 11/22/2024
# author: Yuezhe Li 
# purpose of this code: to generate figure for T-Dxd

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)

using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using Plots, Statistics, DataFrames, CSV, DataFramesMeta

# read observed 
pk_doi =  CSV.read(@projectroot("data/doi-2017.csv"),DataFrame);

# simulation 
include(@projectroot("model/jones_homo.jl") );
include(@projectroot("model/param-pk.jl")); 
include(@projectroot("model/init-pk.jl")); 
include(@projectroot("script/helper.jl")); 

p_dxd = deepcopy(p_base);
# update pk parameters (same as T-DM1)
p_dxd.PS_Score = 6.
p_dxd.init_sR = 0.004
# T-Dxd unique parameters; https://pubmed.ncbi.nlm.nih.gov/37787918/ 
p_dxd.DAR = 8.
p_dxd.k_out = 32.32
p_dxd.k_PL_ex = 0.82

tspan = (-0.01, DayToHour*63);      # [hr]
AddDose = [0., 22., 44.] * DayToHour  # [hr]

iv_dose = [0.8, 1.6, 3.2, 5.4, 6.4, 8.0]; # [mg/kg]

df_tdxd = DataFrame(time = [], pl_lung = [], pl_liver = [], pl_heart = [], pl_muscle = [], pl_skin = [], pl_adipose = [], pl_bone = [], pl_brain = [], 
                  pl_kidney = [], pl_si = [], pl_li = [], pl_pancreas = [], pl_thymus = [], pl_spleen = [], adcplasma = [], dose = []);

for init_dose in iv_dose
    tmpdf = InfusionDoses(init_dose, AddDose, p_dxd, infusion_time = 1);
    tmpdf.dose .= string(init_dose) * "mg/kg"
    df_tdxd = vcat(df_tdxd, tmpdf)
end

# PK
p_pk_tdxd = plot(legend = :top, ylims = [1E-1, 1E4], yaxis = :log, size=(350,350), legendcolumns=3, palette = :Set2_6, dpi = 1000);
plot!(df_tdxd.time/DayToHour, df_tdxd.adcplasma*MW/1E3, group = df_tdxd.dose, linewidth = 2); 
scatter!(pk_doi.time_day, pk_doi.T_DXd_ugperml, group = pk_doi.Dose, ma = 0.6, label = false, linewidth=0);
xlabel!("Time (day)"); xlims!(0, 63); xticks!(0:7:63); ylabel!("plasma T-Dxd concentration (ug/mL)"); 
display(p_pk_tdxd)

savefig(p_pk_tdxd, @projectroot("deliv/figure/pk/t-dxd-homo.png"));

# PAYLOAD 

AddDose_q3w = [0., 21., 42.] * DayToHour  # [hr]

df_tdxd = DataFrame(time = [], pl_lung = [], pl_liver = [], pl_heart = [], pl_muscle = [], pl_skin = [], pl_adipose = [], pl_bone = [], pl_brain = [], 
                  pl_kidney = [], pl_si = [], pl_li = [], pl_pancreas = [], pl_thymus = [], pl_spleen = [], adcplasma = [], dose = []);

for init_dose in iv_dose
    tmpdf = InfusionDoses(init_dose, AddDose_q3w, p_dxd, infusion_time = 1);
    tmpdf.dose .= string(init_dose) * "mg/kg"
    df_tdxd = vcat(df_tdxd, tmpdf)
end

CSV.write(@projectroot("data/sims/homo-tdxd.csv"), df_tdxd)


plung = plot(df_tdxd.time, df_tdxd.pl_lung, group = df_tdxd.dose, title = "lung endothelial payload, tdxd")
pliver = plot(df_tdxd.time, df_tdxd.pl_liver, group = df_tdxd.dose, title = "liver endothelial payload, T-Dxd", xticks = (0:7:63))

savefig(pliver, @projectroot("deliv/figure/payload/t-dxd-q3w-liver.png"));
