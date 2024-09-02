"""
    power(a,b)

    Raise `a` to the `b` power 
"""
function power(a,b)
    a^b
  end 

"""
    ChainsCustomIndex(c::Chains, indices_to_keep::BitMatrix)

Reduce Chains object to only wanted indices. 

Function created by Damon Bayer. 
"""
function ChainsCustomIndex(c::Chains, indices_to_keep::BitMatrix)
    min_length = minimum(mapslices(sum, indices_to_keep, dims = 1))
  v = c.value
  new_v = copy(v.data)
  new_v_filtered = cat([new_v[indices_to_keep[:, i], :, i][1:min_length, :] for i in 1:size(v, 3)]..., dims = 3)
  aa = AxisArray(new_v_filtered; iter = v.axes[1].val[1:min_length], var = v.axes[2].val, chain = v.axes[3].val)

  Chains(aa, c.logevidence, c.name_map, c.info)
end

# Series of functions for creating correctly scaled parameter draws. 
# code snippet shared by @torfjelde
# https://gist.github.com/torfjelde/37be5a672d29e473983b8e82b45c2e41
generate_names(val) = generate_names("", val)
generate_names(vn_str::String, val::Real) = [vn_str;]
function generate_names(vn_str::String, val::NamedTuple)
    return map(keys(val)) do k
        generate_names("$(vn_str)$(k)", val[k])
    end
end
function generate_names(vn_str::String, val::AbstractArray{<:Real})
    results = String[]
    for idx in CartesianIndices(val)
        s = join(idx.I, ",")
        push!(results, "$vn_str[$s]")
    end
    return results
end

function generate_names(vn_str::String, val::AbstractArray{<:AbstractArray})
    results = String[]
    for idx in CartesianIndices(val)
        s1 = join(idx.I, ",")
        inner_results = map(f("", val[idx])) do s2
            "$vn_str[$s1]$s2"
        end
        append!(results, inner_results)
    end
    return results
end

function generate_names(vn_str::String, val::Vector{Any})
    # Added by Christian Bernal Zelaya.
    results = String[]
    for idx in 1:length(val)
        push!(results, "$(vn_str)[$idx]")
    end
    return results
end

flatten(val::Real) = [val;]
function flatten(val::AbstractArray{<:Real})
    return mapreduce(vcat, CartesianIndices(val)) do i
        val[i]
    end
end
function flatten(val::AbstractArray{<:AbstractArray})
    return mapreduce(vcat, CartesianIndices(val)) do i
        flatten(val[i])
    end
end
function flatten(val::Vector{Any})
    # Added by Christian Bernal Zelaya.
    results = []
    for item in val
        append!(results, flatten(item))
    end
    return results
end


function vectup2chainargs(ts::AbstractVector{<:NamedTuple})
    ks = keys(first(ts))
    vns = mapreduce(vcat, ks) do k
        generate_names(string(k), first(ts)[k])
    end
    vals = map(eachindex(ts)) do i
        mapreduce(vcat, ks) do k
            flatten(ts[i][k])
        end
    end
    arr_tmp = reduce(hcat, vals)'
    arr = reshape(arr_tmp, (size(arr_tmp)..., 1)) # treat as 1 chain
    return Array(arr), vns
end

function vectup2chainargs(ts::AbstractMatrix{<:NamedTuple})
    num_samples, num_chains = size(ts)
    res = map(1:num_chains) do chain_idx
        vectup2chainargs(ts[:, chain_idx])
    end

    vals = getindex.(res, 1)
    vns = getindex.(res, 2)

    # Verify that the variable names are indeed the same
    vns_union = reduce(union, vns)
    @assert all(isempty.(setdiff.(vns, Ref(vns_union)))) "variable names differ between chains"

    arr = cat(vals...; dims = 3)

    return arr, first(vns)
end

function MCMCChains.Chains(ts::AbstractArray{<:NamedTuple})
    return MCMCChains.Chains(vectup2chainargs(ts)...)
end


"""
    save_plots_to_docs(plot, filename; format = "png")

Saves plots to docs/plots directory.

"""
function save_plots_to_docs(plot, filename; format = "png")
    doc_loc = "plots"
    if !isdir(doc_loc)
        mkdir(doc_loc)
    end

    file_target_path = joinpath(doc_loc, "$filename.$format")
    savefig(plot, file_target_path)
    println("Plot saved to $file_target_path")

end


"""
    startswith_any(name, patterns)

Checks if the name of time varying paramter starts with any of the patterns.

"""
function startswith_any(name, patterns)
    for pattern in patterns
        if startswith(name, pattern)
            return true
        end
    end
    return false
end


"""
    calculate_quantiles(df, chain, var_prefix, quantiles)

Calculate quantiles for a given chain and variable prefix.  Quantiles can be any user desired quantile.

"""
function calculate_quantiles(df, chain, var_prefix, quantiles)
    df_chain = filter(row -> row.chain == chain, df)
    column_names = names(df_chain)
    var_names = filter(name -> startswith_any(name, [var_prefix]), column_names)
    medians = [median(df_chain[:, var]) for var in var_names]  
    lower_bounds = [quantile(df_chain[:, var], (1 .- quantiles) / 2) for var in var_names]
    upper_bounds = [quantile(df_chain[:, var], 1 .- (1 .- quantiles) / 2) for var in var_names]
    
    
    return medians, lower_bounds, upper_bounds
end


"""
    generate_ribbon_colors(number_of_colors)

Generates a vector with colors for ribbons in plots.

"""
function generate_colors(number_of_colors)
    alpha_values = range(0.1, stop=0.7, length=number_of_colors)
    return [RGBA(colorant"blue", alpha) for alpha in alpha_values]
end


"""
    repeat_last_n_elements(x::Vector{T}, n::Int, w::Int) where T

Modifies a given array so that the last n elements are repeated w times.
    
"""
function repeat_last_n_elements(x::Vector{T}, n::Int, w::Int) where T
    if n == 0
        return x
    else
        n = min(n, length(x))
        last_n_elements = x[end-n+1:end]
        repeated_elements = [elem for elem in last_n_elements for _ in 1:w]
        x_new = vcat(x, repeated_elements)
    
        return x_new
    end
end