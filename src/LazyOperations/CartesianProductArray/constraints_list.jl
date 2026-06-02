"""
    constraints_list(cpa::CartesianProductArray)

Compute a list of constraints of a (polyhedral) Cartesian product of a finite
number of sets.

### Input

- `cpa` -- Cartesian product of a finite number of sets

### Output

A list of constraints.
"""
function constraints_list(cpa::CartesianProductArray)
    return _constraints_list_cartesian_product(cpa)
end

function _constraints_list_cartesian_product(cp::Union{CartesianProduct,CartesianProductArray})
    N = eltype(cp)
    clist = Vector{HalfSpace{N,SparseVector{N,Int}}}()
    n = dim(cp)
    sizehint!(clist, n)
    prev_step = 1
    # create high-dimensional constraints list
    for c_low in cp
        c_low_list = constraints_list(c_low)
        if isempty(c_low_list)
            n_low = dim(c_low)
        else
            n_low = dim(c_low_list[1])
            indices = prev_step:(prev_step + n_low - 1)
        end
        for constr in c_low_list
            new_constr = HalfSpace(sparsevec(indices, constr.a, n), constr.b)
            push!(clist, new_constr)
        end
        prev_step += n_low
    end

    return clist
end
