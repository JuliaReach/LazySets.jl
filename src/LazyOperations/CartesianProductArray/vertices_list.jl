"""
    vertices_list(cpa::CartesianProductArray)

Compute a list of vertices of a (polytopic) Cartesian product of a finite
number of sets.

### Input

- `cpa` -- Cartesian product of a finite number of sets

### Output

A list of vertices.

### Algorithm

We assume that the underlying sets are polytopic.
Then the high-dimensional set of vertices is just the Cartesian product of the
low-dimensional sets of vertices.
"""
@validate function vertices_list(cpa::CartesianProductArray)
    # collect low-dimensional vertices lists
    vlist_low = [vertices_list(X) for X in cpa]

    # create high-dimensional vertices list
    indices_max = [length(vl) for vl in vlist_low]
    m = prod(indices_max)
    N = eltype(cpa)
    vlist = Vector{Vector{N}}(undef, m)
    indices = ones(Int, length(vlist_low))
    v = zeros(N, dim(cpa))
    dim_start_j = 1
    for vl in vlist_low
        v_low = vl[1]
        v[dim_start_j:(dim_start_j + length(v_low) - 1)] = v_low
        dim_start_j += length(v_low)
    end
    i = 1
    j = 1
    # iterate through all index combinations
    while true
        indices[1] = 0
        while indices[1] < indices_max[1]
            indices[1] += 1
            v_low = vlist_low[1][indices[1]]
            v[1:length(v_low)] = v_low
            vlist[i] = copy(v)
            i += 1
        end
        if i > m
            break
        end
        j = 1
        dim_start_j = 1
        while indices[j] == indices_max[j]
            indices[j] = 1
            v_low = vlist_low[j][1]
            v[dim_start_j:(dim_start_j + length(v_low) - 1)] = v_low
            dim_start_j += length(v_low)
            j += 1
        end
        indices[j] += 1
        v_low = vlist_low[j][indices[j]]
        v[dim_start_j:(dim_start_j + length(v_low) - 1)] = v_low
    end

    return vlist
end
