# date: 7/1/2025 
# author: Yuezhe Li 
# purpose of this code: to save outcome 

function ProcessOutcome(sol_dict)
    df = DataFrame(time_hr = [], adc_plasma_uM = [], pl_liver_uM = [], dose_mgkg = [])

    dose_mgkg_list = collect(keys(sol_dict))

    for dose_i in dose_mgkg_list
        tmpdf = DataFrame(
            time_hr = sol_dict[dose_i].t, 
            adc_plasma_uM = [sol_dict[dose_i].u[i].C_EXG_Plasma for i in 1:length(sol_dict[dose_i].t)],  # [uM]
            pl_liver_uM = [sol_dict[dose_i].u[i].end_cyto_payload[2] for i in 1:length(sol_dict[dose_i].t)],  # [uM]
        ); 
        tmpdf.dose_mgkg .= dose_i; # [mg/kg]
        df = vcat(df, tmpdf)
    end
    @rsubset!(df, :time_hr >=0)
    return df 
end