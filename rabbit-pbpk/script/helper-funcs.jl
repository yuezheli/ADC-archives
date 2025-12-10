###

## helper functions

### Tim's functions

## get parameter map

function get_p_map(model, organ, parameter, value)
    param_array = filter(var -> occursin(organ, string(var)) && occursin(parameter, string(var)), parameters(model))
    p_map = []
    for x in param_array
        push!(p_map, Pair(x, value))
    end
    return vcat(p_map...)
end

function get_p_map_all(model, p_dict)
    param_map = []
    organs = keys(p_dict)

    for x in organs
        parval = p_dict[x]
        organ_map = []
        for j in 1:length(parval)
            param = string(keys(parval)[j])
            value = values(parval)[j]
            param_new = get_p_map(model, x, param, value)
            organ_map = push!(organ_map, param_new)
        end
        organ_map_flat = vcat(organ_map...)
        param_map = push!(param_map, organ_map_flat)
    end
    param_map_flat = vcat(param_map...)

    return param_map_flat
end
