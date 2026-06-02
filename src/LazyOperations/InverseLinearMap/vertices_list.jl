"""
    vertices_list(ilm::InverseLinearMap; prune::Bool=true)

Return the list of vertices of a (polyhedral) inverse linear map.

### Input

- `ilm`   -- inverse linear map
- `prune` -- (optional, default: `true`) if `true`, remove redundant vertices

### Output

A list of vertices.

### Algorithm

We assume that the underlying set `X` is polyhedral.
Then the result is just the inverse linear map applied to the vertices of `X`.
"""
@validate function vertices_list(ilm::InverseLinearMap; prune::Bool=true)
    # collect vertices list of wrapped set
    vlist_X = vertices_list(ilm.X)

    # create resulting vertices list
    vlist = Vector{eltype(vlist_X)}(undef, length(vlist_X))
    @inbounds for (i, vi) in enumerate(vlist_X)
        vlist[i] = ilm.M \ vi
    end

    return prune ? convex_hull(vlist) : vlist
end
