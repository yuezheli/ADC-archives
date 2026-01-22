# date: 7/1/2025 
# author: Yuezhe Li 
# purpose of this code: to analyze biodistribution of an ADC (use T-DM1 as example)

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)
using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using Parameters: @unpack
using DataFrames
using Plots

include(@projectroot("Constants.jl"))
include(@projectroot("model/jones_homo.jl") )
include(@projectroot("model/param-pk.jl"))
include(@projectroot("model/init-pk.jl")) 
include(@projectroot("script/helper-infusion-dosing.jl"))
include(@projectroot("script/helper-human-organ-volumes.jl"))
include(@projectroot("script/helper-tissue-mass.jl"))

# set dose and time for simulation 
tspan = (0., hr_per_day*84);      # [hr]
Dose = 3.6  # [mg/kg]

dose_umol = Dose*BW*1E-3/MW_IGG*1E6; # [umol]

# set TMDD in plasma to 0 
p_notmdd = deepcopy(p_base); 
p_notmdd.init_sR = 0; 

# simulation (IV bolus dosing)
u0_tmp = jones_init(Dose*1E3, p_notmdd, BW, V_Plasma);
sol_tmp = solve(ODEProblem(jonesODEs_homo_tumor!, u0_tmp, tspan, p_notmdd), saveat = 0.25, alg = QNDF(), reltol = 1E-12);

# ADC ended up in a tissue
function total_tissue_adc(sol; N_Organs = 15)
    # total ADC in tissue (umol)
    tm_adc = HumanTissueMass(sol, N_Organs = N_Organs); 

    # degraded ADC 
    deg_adc = zeros(length(sol.t), N_Organs);
    for i in 1:N_Organs
        deg_adc[:, i] = [sol.u[j].DEGprotein.deg_EXG[i] for j in 1:length(sol.t)]
    end

    # compute total adc 
    total_adc_umol = tm_adc .+ deg_adc

    return total_adc_umol
end

biodist_adc = total_tissue_adc(sol_tmp, N_Organs = 16); 

Organ_Names = ["Lung", "Liver", "Heart", "Muscle", "Skin", "Adipose", "Bone", "Brain", "Kidney", 
                "Small Intestin", "Large Intestin",  "Pancreas", "Thymus", "Spleen", "Other", "Bone Marrow"];

dist_adc_d63 = DataFrame(
    Organ = Organ_Names, 
    frac =  biodist_adc[end, :] / dose_umol,
)

sort!(dist_adc_d63, [:frac], rev = true)

p_total_mass_tissue = bar(dist_adc_d63.Organ, dist_adc_d63.frac * 100, label = false,
    xticks =:all, xrotation = 45, xlabel ="", ylabel = "% of the total drug", xtickfontsize=10, 
    ylims = (5E-3, 15), yaxis = :log, yticks = ([5E-3, 1E-2, 1E-1, 1, 5, 10, 15], ["0.005", "0.01", "0.1", "1", "5","10", "15"]),
    size = [400, 400], dpi = 300);

savefig(p_total_mass_tissue, @projectroot("deliv/figure/biodistribution/t-dm1.png"));

# Cmax of ADCs in tissue endothelial cells 
function tissue_adc_cmax(sol; N_Organs = 15)

    # total ADC in tissue (umol)
    tm_adc = HumanTissueMass(sol, N_Organs = N_Organs); 

    # organ volume 
    V_Organs = TotalHumanOrganVolume(N_Organs = N_Organs); 

    # Concentration 
    tissue_adc_conc = zeros(size(tm_adc))
    for i in 1:size(tm_adc)[1]
        tissue_adc_conc[i, :] = tm_adc[i, :] ./ V_Organs
    end
    
    # Cmax 
    Cmax = []
    for i in 1:size(tm_adc)[2]
        tmp_conc = tissue_adc_conc[:, i]
        tmp_cmax = maximum(tmp_conc)
        append!(Cmax, tmp_cmax)
    end

    return Cmax
end

Cmax_ADCs = tissue_adc_cmax(sol_tmp, N_Organs = 16);

cmax_adc_d63 = DataFrame(
    Organ = Organ_Names, 
    cmax_nM = Cmax_ADCs * 1E3
); 

sort!(cmax_adc_d63, [:cmax_nM], rev = true)

p_cmax_tissue = bar(cmax_adc_d63.Organ, cmax_adc_d63.cmax_nM, label = false,
    xticks =:all, xrotation = 45, xlabel ="", ylabel = "Maximum ADC concentration(nM)", xtickfontsize=10, 
    ylims = (5, 150), yaxis = :log, yticks = ([5, 10, 50, 100, 150], ["5", "10", "50","100", "150"]),
    size = [400, 400], dpi = 300);

savefig(p_cmax_tissue, @projectroot("deliv/figure/biodistribution/t-dm1-cmax.png"));

