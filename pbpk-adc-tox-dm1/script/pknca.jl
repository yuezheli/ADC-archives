# author: Yuezhe Li 
# date: 4/3/24
# purpose of this script: functions for noncompartmental PK analysis

using Pkg; Pkg.activate("");

# area under the curve 
function auc_sparse(time, conc)
    tmp = 0; 
    if length(time) == length(conc)
        for i = 1:(length(time)-1)
            t1 = time[i]
            t2 = time[i+1]
            c1 = conc[i]
            c2 = conc[i+1]
            tmparea = (c1 + c2)*(t2-t1)/2
            tmp += tmparea
        end
        return tmp
    else
        println("wrong input")
        return nothing; 
    end
end

# Area under the Moment Curve
function aumc_sparse(time, conc)
    # time, conc should all be vectors
    tmp = 0; 
    if length(time) == length(conc)
        for i = 1:(length(time)-1)
            t1 = time[i]
            t2 = time[i+1]
            c1 = conc[i]
            c2 = conc[i+1]
            tmparea = (t1*c1 + t2*c2)*(t2-t1)/2
            tmp += tmparea
        end
        return tmp
    else
        println("wrong input")
        return nothing; 
    end
end

# mean residence time; coded following https://www.certara.com/knowledge-base/mean-residence-time-mrt-understanding-how-long-drug-molecules-stay-in-the-body/
# note the function is only good for bolus IV dosing 
function mrt_sparse(time, conc)
    return aumc_sparse(time, conc)/auc_sparse(time, conc)
end

# steady state volume of distribution; coded based on https://www.rdocumentation.org/packages/PKNCA/versions/0.9.1/topics/pk.calc.vss
function Vss_sparse(dose, time, conc)
    cl = dose/auc_sparse(time, conc);
    vol_ss = cl * mrt_sparse(time, conc);
    return vol_ss
end

# calculate half life 
function half_life_sparse(t, conc)
    i = 1;
    try 
        indx = findall(x->x==maximum(conc), conc);
        i = indx[end]; 
    catch 
        i = 1; 
    end
    while i < length(t)
        if conc[i] < maximum(conc)/2
            break
        else
            i += 1
        end
    end
    return t[i]
end
