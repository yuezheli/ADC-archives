# date: 7/2/2025 
# author: Yuezhe Li 
# purpose of this: to calculate t1/2 

function find_closest_index(list, target)
    abs_diffs = abs.(list .- target)
    return argmin(abs_diffs)
end

function plasma_half_life(time, plasmaconc; ivbolus = true, todisgard = 2)
    if ivbolus
        cmax = maximum(plasmaconc[todisgard:end])
        chalf = cmax/2
        thalf_index = find_closest_index(plasmaconc, chalf)
        thalf = time[thalf_index]
    else 
        cmax = maximum(plasmaconc)
        chalf = cmax/2
        chalf = cmax/2
        thalf_index = find_closest_index(plasmaconc, chalf)
        thalf = time[thalf_index]
    end
    return thalf
end
