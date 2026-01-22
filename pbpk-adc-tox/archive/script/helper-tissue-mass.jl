# date: 7/1/2025 
# author: Yuezhe Li 
# purpose of this code: to calculate the amount of ADC (umol) distributed into organ (a sum of vascular, membrane, cellular)

function HumanTissueMass(sol; N_Organs = 16)
    V_V, V_VM, V_E7, V_E6a, V_E7b, V_ISM, V_IntS = ReturnHumanOrganVolume(N_Organs = N_Organs)

    # mAb Location Indices
    V = 1
    VM = 2
    E7 = 3
    E6a = 4
    E7b = 5
    ISM = 6
    IntS = 7
    bound_VM = 8
    bound_E7 = 9
    bound_E6a = 10
    bound_E7b = 11
    bound_ISM = 12
    bound2_VM = 13
    bound2_E7 = 14
    bound2_E6a = 15
    bound2_E7b = 16
    bound2_ISM = 17
    bound_VM_mem = 18
    bound_ISM_mem = 19

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
 
    ## sum up mAb mass in each organ
    tmp2 = zeros(length(sol.t), N_Organs);
    for i = 1:length(sol.t)
        for j = 1:N_Organs
            tmp2[i,j] = sum(tmp[i,j,:])
        end
    end

    return tmp2;
end