# date: 1/27/2025
# author: Yuezhe Li 
# purpose of this code: to generate simulation of Cantuzumab mertansine PK 

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)
using DataFrames, CSV, DataFramesMeta

# read observed; https://pubmed.ncbi.nlm.nih.gov/12525512/; The dose-limiting toxicity (DLT) of cantuzumab mertansine was found to be reversible elevations of hepatic transaminases
pk_tolcher =  CSV.read("data/Tolcher-2003.csv",DataFrame);

using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using Plots

# simulation 
include(@projectroot("model/jones_homo.jl") );
include(@projectroot("model/param-pk.jl")); 
include(@projectroot("model/init-pk.jl")); 
include(@projectroot("script/helper.jl"));

# DAR value not changed based on https://pubmed.ncbi.nlm.nih.gov/18301896/
p_cm = deepcopy(p_base);
p_cm.PS_Score = -2
p_cm.PS_kd = 0.0001
p_cm.init_sR = 0.001  

tspan = (-0.01, DayToHour*84);      # [hr]
AddDose_q3w = [0., 21., 42., 63] * DayToHour  # [hr]

q3w_dose = [22, 44, 88, 132, 176, 235, 295]; # [mg/m2]

df_cm_q3w = DataFrame(time = [], pl_lung = [], pl_liver = [], pl_heart = [], pl_muscle = [], pl_skin = [], pl_adipose = [], pl_bone = [], pl_brain = [], 
                  pl_kidney = [], pl_si = [], pl_li = [], pl_pancreas = [], pl_thymus = [], pl_spleen = [], adcplasma = [], dose = []);

for init_dose in q3w_dose
    tmpdf = InfusionDoses(init_dose*HT*HT/BW, AddDose_q3w, p_cm, infusion_time = 0.5);
    tmpdf.dose .= string(init_dose) * " mg/m2 Q3W"
    df_cm_q3w = vcat(df_cm_q3w, tmpdf)
end

df_cm_235 = @rsubset(df_cm_q3w, :dose == "235 mg/m2 Q3W");

# PK dynamics 
p_pk_cm = plot(legend = :topright, ylims = [1E2, 1E6], xlims = [0, 500], yaxis = :log, size=(350,350), dpi = 1000);
plot!(df_cm_235.time, df_cm_235.adcplasma*MW, linewidth = 2, alpha = 0.8, color = "red", label = "235mg/m² Q3W"); 
scatter!(pk_tolcher.time_hour, pk_tolcher.conc_ug_L, ma = 0.6, label = false, linewidth=0);
xlabel!("Time (hour)"); ylabel!("plasma cantuzumab mertansine (ug/L)", guidefontsize = 10); 
display(p_pk_cm)

savefig(p_pk_cm, @projectroot("deliv/figure/pk/cantuzumab-mertansine-homo.png"));

# payload 
pliver = plot(df_cm_q3w.time/DayToHour, df_cm_q3w.pl_liver, group = df_cm_q3w.dose, title = "liver endothelial payload, cantuzumab mertansine", ylabel = "DM1 concentration (uM)", titlefontsize = 8, xticks = (0:7:84))

savefig(pliver, @projectroot("deliv/figure/payload/cantuzumab-mertansine-q3w-liver.png"));

# additional dosing scheme from https://aacrjournals.org/clincancerres/article/10/13/4363/94515/A-Phase-I-Study-of-Cantuzumab-Mertansine
# dosing QW; DLT is liver tox, at 115 mg/m2 QW

qw_dose = [40, 60, 80, 96, 115, 138]; # [mg/m2]
AddDose_qw = [0., 7., 14., 21., 28., 35., 42., 49., 56., 63., 70., 77.] * DayToHour  # [hr]

df_cm_qw = DataFrame(time = [], pl_lung = [], pl_liver = [], pl_heart = [], pl_muscle = [], pl_skin = [], pl_adipose = [], pl_bone = [], pl_brain = [], 
                  pl_kidney = [], pl_si = [], pl_li = [], pl_pancreas = [], pl_thymus = [], pl_spleen = [], adcplasma = [], dose = []);

for init_dose in qw_dose
    tmpdf = InfusionDoses(init_dose*HT*HT/BW, AddDose_qw, p_cm, infusion_time = 0.5);
    tmpdf.dose .= string(init_dose) * " mg/m2 QW"
    df_cm_qw = vcat(df_cm_qw, tmpdf)
end

# add another dosing scheme from Rodin 2008 
# https://link.springer.com/article/10.1007/s00280-007-0672-8

qw34_dose = [30, 45, 60]; # [mg/m2]
AddDose_qw34 = [0, 2, 4, 7, 9, 11, 14, 16, 18, 
                28, 30, 32, 35, 37, 39, 42, 44, 46, 
                56, 58, 60, 63, 65, 67, 70, 72, 74] * DayToHour  # [hr]

df_cm_qw34 = DataFrame(time = [], pl_lung = [], pl_liver = [], pl_heart = [], pl_muscle = [], pl_skin = [], pl_adipose = [], pl_bone = [], pl_brain = [], 
                  pl_kidney = [], pl_si = [], pl_li = [], pl_pancreas = [], pl_thymus = [], pl_spleen = [], adcplasma = [], dose = []);

for init_dose in qw34_dose
    tmpdf = InfusionDoses(init_dose*HT*HT/BW, AddDose_qw34, p_cm, infusion_time = 0.5);
    tmpdf.dose .= string(init_dose) * " mg/m2 QOD"
    df_cm_qw34 = vcat(df_cm_qw34, tmpdf)
end

# add another dose that matched 3.6 mg/kg 

macthed_df = InfusionDoses(3.6, AddDose_q3w, p_cm, infusion_time = 0.5);
macthed_df.dose .= "3.6 mg/kg Q3W";

df_cm = vcat(df_cm_q3w, df_cm_qw, df_cm_qw34, macthed_df);

CSV.write(@projectroot("data/sims/homo-cantuzumab-mertansine.csv"), df_cm)
