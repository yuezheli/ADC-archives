# author: Yuezhe LI
# date: 

using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using Statistics, DataFrames, CSV, DataFramesMeta
using Parameters: @unpack

include("params.jl") # load this for number of organs in the model 

# create initial value for human status
function solid_tumor_init(N_Organs, TotalCells, RposFrac, Dose_in_ugkg, Rcopy, init_sR, k_endo, k_rec, isnodal = false)

    InitRpos = TotalCells * RposFrac
    InitRneg = TotalCells * (1-RposFrac)

    MW_EDG = 150000.0; # [Da]
    FcRn_Conc = 44.093; # [uM]  
    EDG_mg_ml = 10.0;  # [mg/mL]
    C_EDG_Plasma0 =  EDG_mg_ml/MW_EDG*1E6; # [uM]
    C_EDG_LN0 = 1E-18*1E6; # [uM]
    C_EDG0 =  1E-18*1E6 * ones(N_Organs,17); # [uM]
    C_EXG_LN0 = 1E-18*1E6; # [uM]
    C_EXG0 =  1E-18*1E6 *ones(N_Organs,19); # [uM]
    C_FcRn_E6a0 =  ones(N_Organs)*FcRn_Conc; # [uM]
    C_FcRn_E7b0 =  ones(N_Organs)*FcRn_Conc; # [uM] 
    C_FcRn_E70 =  ones(N_Organs)*FcRn_Conc; # [uM]      
    C_FcRn_ISM0 =  ones(N_Organs)*FcRn_Conc*1E-4; # [uM] 
    C_FcRn_VM0 =  ones(N_Organs)*FcRn_Conc*1E-4; # [uM] 
    C_EXG_Plasma0 =  Dose_in_ugkg*BW/(V_Plasma)/MW_EDG; # [uM]

    init_DEGprotein = ComponentArray(deg_EXG = zeros(N_Organs), deg_EDG = zeros(N_Organs), deg_TumorCellular = 0.0, deg_plasma = 0.0)

    if isnodal
        # this is for nodal tumor
        init_tumor = ComponentArray(
        Nc_1 = InitRpos, Nc_2 = 0., Nc_3 = 0., Nc_4 = 0., Nc_1_neg = InitRneg, Nc_2_neg = 0., Nc_3_neg = 0., Nc_4_neg = 0.,
        R_s = Rcopy/N_av*1e6, R_e = k_endo/k_rec*Rcopy/N_av*1e6, AR_s = 0., AR_e = 0., P_c = 0., P_m = 0., P_neg_c = 0.);
    else
        # these are solid tumor
        if N_Organs == 15
            # this is when tumor was modeled following Shah et al., 2012
            init_tumor = ComponentArray(
            Nc_1 = InitRpos, Nc_2 = 0., Nc_3 = 0., Nc_4 = 0., 
            Nc_1_neg = InitRneg, Nc_2_neg = 0., Nc_3_neg = 0., Nc_4_neg = 0.,
            R_s = Rcopy/N_av*1e6, R_e = k_endo/k_rec*Rcopy/N_av*1e6, AR_s = 0., AR_e = 0., P_c = 0., P_m = 0., P_neg_c = 0., 
            A_m = 0.
            );
        else
            # this is when tumor was modeled as another organ 
            init_tumor = ComponentArray(
            Nc_1 = InitRpos, Nc_2 = 0., Nc_3 = 0., Nc_4 = 0., Nc_1_neg = InitRneg, Nc_2_neg = 0., Nc_3_neg = 0., Nc_4_neg = 0.,
            R_s = Rcopy/N_av*1e6, R_e = k_endo/k_rec*Rcopy/N_av*1e6, AR_s = 0., AR_e = 0., P_c = 0., P_m = 0., P_neg_c = 0.);
        end
    end

    u0 = ComponentArray(C_EXG_Plasma = C_EXG_Plasma0, 
                            C_sR_plasma = init_sR, C_sR_EXG = 0., 
                            C_EDG_Plasma = C_EDG_Plasma0, C_EDG_LN = C_EDG_LN0, C_EDG = C_EDG0, C_EXG_LN = C_EXG_LN0,
                            C_EXG = C_EXG0, C_FcRn_E6a = C_FcRn_E6a0, C_FcRn_E7 = C_FcRn_E70, C_FcRn_E7b = C_FcRn_E7b0, 
                            C_FcRn_ISM = C_FcRn_ISM0, C_FcRn_VM = C_FcRn_VM0,tumor = init_tumor, DEGprotein = init_DEGprotein, 
                            end_endo_payload = zeros(N_Organs), end_cyto_payload = zeros(N_Organs), ints_payload = zeros(N_Organs)); 

    return u0, Dose_in_ugkg*BW/MW_EDG ;
end


# this function aimed to add up all ADC in each organ (in plasma, in interstitium, and all parts of endothelial cells); default unit in umol.
function TissueMass(sol, V_V, V_VM, V_E7, V_E6a, V_E7b, V_ISM, V_IntS)
    # this function calculates the total mass of ADC inside tissue (vascular, membrane, cellular)
    V = 1; VM = 2; E7 = 3; E6a = 4; E7b = 5; ISM = 6; IntS = 7; 
    bound_VM = 8; bound_E7 = 9; bound_E6a = 10; bound_E7b = 11; bound_ISM = 12
    bound2_VM = 13; bound2_E7 = 14; bound2_E6a = 15; bound2_E7b = 16; bound2_ISM = 17
    bound_VM_mem = 18; bound_ISM_mem = 19

    N_Organs = length(V_V)

    ## extract mAb mass from different parts of each organ
    tmp = zeros(length(sol.t), N_Organs, 19);

    tmp[:, :, V] = mapreduce(permutedims, vcat, [sol.u[i].C_EXG[1:N_Organs,V].* V_V[1:N_Organs]  for i in 1:length(sol.t)] )
    tmp[:, :, VM] = mapreduce(permutedims, vcat, [sol.u[i].C_EXG[1:N_Organs,VM] .* V_VM[1:N_Organs]  for i in 1:length(sol.t)] )
    tmp[:, :, E7] = mapreduce(permutedims, vcat, [sol.u[i].C_EXG[1:N_Organs,E7] .* V_E7[1:N_Organs]  for i in 1:length(sol.t)] )
    tmp[:, :, E6a] = mapreduce(permutedims, vcat, [sol.u[i].C_EXG[1:N_Organs,E6a] .* V_E6a[1:N_Organs]  for i in 1:length(sol.t)] )
    tmp[:, :, E7b] = mapreduce(permutedims, vcat, [sol.u[i].C_EXG[1:N_Organs,E7b] .* V_E7b[1:N_Organs]  for i in 1:length(sol.t)] )
    tmp[:, :, ISM] = mapreduce(permutedims, vcat, [sol.u[i].C_EXG[1:N_Organs,ISM] .* V_ISM[1:N_Organs]  for i in 1:length(sol.t)] )
    tmp[:, :, IntS] = mapreduce(permutedims, vcat, [sol.u[i].C_EXG[1:N_Organs,IntS] .* V_IntS[1:N_Organs]  for i in 1:length(sol.t)] )

    tmp[:, :, bound_VM] = mapreduce(permutedims, vcat, [sol.u[i].C_EXG[1:N_Organs,bound_VM] .* V_VM[1:N_Organs]  for i in 1:length(sol.t)] )
    tmp[:, :, bound_E7] = mapreduce(permutedims, vcat, [sol.u[i].C_EXG[1:N_Organs,bound_E7] .* V_E7[1:N_Organs]  for i in 1:length(sol.t)] )
    tmp[:, :, bound_E6a] = mapreduce(permutedims, vcat, [sol.u[i].C_EXG[1:N_Organs,bound_E6a] .* V_E6a[1:N_Organs]  for i in 1:length(sol.t)] )
    tmp[:, :, bound_E7b] = mapreduce(permutedims, vcat, [sol.u[i].C_EXG[1:N_Organs, bound_E7b] .* V_E7b[1:N_Organs]  for i in 1:length(sol.t)] )
    tmp[:, :, bound_ISM] = mapreduce(permutedims, vcat, [sol.u[i].C_EXG[1:N_Organs,bound_ISM] .* V_ISM[1:N_Organs]  for i in 1:length(sol.t)] )

    tmp[:, :, bound2_VM] = mapreduce(permutedims, vcat, [sol.u[i].C_EXG[1:N_Organs,bound2_VM] .* V_VM[1:N_Organs]  for i in 1:length(sol.t)] )
    tmp[:, :, bound2_E7] = mapreduce(permutedims, vcat, [sol.u[i].C_EXG[1:N_Organs,bound2_E7] .* V_E7[1:N_Organs]  for i in 1:length(sol.t)] )
    tmp[:, :, bound2_E6a] = mapreduce(permutedims, vcat, [sol.u[i].C_EXG[1:N_Organs,bound2_E6a] .* V_E6a[1:N_Organs]  for i in 1:length(sol.t)] )
    tmp[:, :, bound2_E7b] = mapreduce(permutedims, vcat, [sol.u[i].C_EXG[1:N_Organs, bound2_E7b] .* V_E7b[1:N_Organs]  for i in 1:length(sol.t)] )
    tmp[:, :, bound2_ISM] = mapreduce(permutedims, vcat, [sol.u[i].C_EXG[1:N_Organs,bound2_ISM] .* V_ISM[1:N_Organs]  for i in 1:length(sol.t)] )

    tmp[:, :, bound_VM_mem] = mapreduce(permutedims, vcat, [sol.u[i].C_EXG[1:N_Organs, bound_VM_mem] .* V_VM[1:N_Organs]  for i in 1:length(sol.t)] )
    tmp[:, :, bound_ISM_mem] = mapreduce(permutedims, vcat, [sol.u[i].C_EXG[1:N_Organs, bound_ISM_mem] .* V_ISM[1:N_Organs]  for i in 1:length(sol.t)] )

    ## calculate ADC inside tumor cells
    TumorCellAb = 
    try
        [ (sol.u[i].tumor.AR_s+sol.u[i].tumor.AR_e+sol.u[i].tumor.A_e+sol.u[i].tumor.AR_l+sol.u[i].tumor.A_l).*(sol.u[i].tumor.Nc_1+sol.u[i].tumor.Nc_2+sol.u[i].tumor.Nc_3+sol.u[i].tumor.Nc_4) for i in 1:length(sol.t)];
    catch 
        # incorporate to deal with a simplied tumor model
        [ (sol.u[i].tumor.AR_s + sol.u[i].tumor.AR_e).*(sol.u[i].tumor.Nc_1) for i in 1:length(sol.t)] ;
    end

    ## sum up mAb mass in each organ
    tmp2 = zeros(length(sol.t), N_Organs);
    for i = 1:length(sol.t)
        for j = 1:N_Organs
            tmp2[i,j] = sum(tmp[i,j,:])
        end
        # tmp2[i,Tumor] += TumorCellAb[i]
    end

    return tmp2;
end

# this function tracks ADCs beingd degraded by organ endothelial cells during the distribution
# default unit: umol
function EpithelialDegradedProtein(sol)
    
    N_Organs = length(V_V)

    tmp = zeros(length(sol.t), N_Organs);

    for Organ in 1:N_Organs
        tmp[:, Organ] = [sol.u[i].DEGprotein.deg_EXG[Organ] for i in 1:length(sol.t)];
    end

    return tmp;
end

# this function calculate average tissue ADC concentration (assuming the tissue of interest was ground up)
function TissueConc(sol)
    totalAb = TissueMass(sol, V_V, V_VM, V_E7, V_E6a, V_E7b, V_ISM, V_IntS); 
    OrganConc = [totalAb[i,:] ./ V_Organ for i in 1:size(totalAb)[1]];
    OrganConc2 = mapreduce(permutedims, vcat, OrganConc);       # convert the previous list to a 2d table
    tumorcells = [ (sol.u[i].tumor.Nc_1+sol.u[i].tumor.Nc_2+sol.u[i].tumor.Nc_3+sol.u[i].tumor.Nc_4+
            sol.u[i].tumor.Nc_1_neg+sol.u[i].tumor.Nc_2_neg+sol.u[i].tumor.Nc_3_neg+sol.u[i].tumor.Nc_4_neg) for i in 1:length(sol.t)];
    return OrganConc2;
end


# extract ADC concentration in tissue interstitium 

function InterstitialADC(sol)
    IntS = 7
    lung_ints =     [sol.u[i].C_EXG[1, IntS] for i in 1:length(sol.t)]; # [uM]
    liver_ints =    [sol.u[i].C_EXG[2, IntS] for i in 1:length(sol.t)]; # [uM]
    heart_ints =    [sol.u[i].C_EXG[3, IntS] for i in 1:length(sol.t)]; # [uM]
    muscle_ints =   [sol.u[i].C_EXG[4, IntS] for i in 1:length(sol.t)]; # [uM]
    skin_ints =     [sol.u[i].C_EXG[5, IntS] for i in 1:length(sol.t)]; # [uM]
    adipose_ints =  [sol.u[i].C_EXG[6, IntS] for i in 1:length(sol.t)]; # [uM]
    bone_ints =     [sol.u[i].C_EXG[7, IntS] for i in 1:length(sol.t)]; # [uM]
    brain_ints =    [sol.u[i].C_EXG[8, IntS] for i in 1:length(sol.t)]; # [uM]
    kidney_ints =   [sol.u[i].C_EXG[9, IntS] for i in 1:length(sol.t)]; # [uM]
    si_ints =       [sol.u[i].C_EXG[10, IntS] for i in 1:length(sol.t)]; # [uM]
    li_ints =       [sol.u[i].C_EXG[11, IntS] for i in 1:length(sol.t)]; # [uM]
    pancreas_ints = [sol.u[i].C_EXG[12, IntS] for i in 1:length(sol.t)]; # [uM]
    thymus_ints =   [sol.u[i].C_EXG[13, IntS] for i in 1:length(sol.t)]; # [uM]
    spleen_ints =   [sol.u[i].C_EXG[14, IntS] for i in 1:length(sol.t)]; # [uM]
    df = DataFrame(
        lung_ints = lung_ints, liver_ints = liver_ints, heart_ints = heart_ints, muscle_ints = muscle_ints, 
        skin_ints = skin_ints, adipose_ints = adipose_ints, bone_ints = bone_ints, brain_ints = brain_ints, 
        kidney_ints = kidney_ints, si_ints = si_ints, li_ints = li_ints, pancreas_ints = pancreas_ints, 
        thymus_ints = thymus_ints, spleen_ints = spleen_ints
    );
    return df; 
end


