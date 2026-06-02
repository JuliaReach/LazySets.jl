"""
    vertices_list(lm::LinearMap; prune::Bool=true)

Return the list of vertices of a (polytopic) linear map.

### Input

- `lm`    -- linear map
- `prune` -- (optional, default: `true`) if `true`, we remove redundant vertices

### Output

A list of vertices.

### Algorithm

We assume that the underlying set `X` is polytopic and compute the vertices of
`X`. The result is just the linear map applied to each vertex.
"""
@validate function vertices_list(lm::LinearMap; prune::Bool=true)
    # apply the linear map to each vertex
    vlist = broadcast(x -> lm.M * x, vertices(lm.X))

    return prune ? convex_hull(vlist) : vlist
end
