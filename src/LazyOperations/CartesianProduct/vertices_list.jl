"""
    vertices_list(cp::CartesianProduct)

Return the list of vertices of a (polytopic) Cartesian product.

### Input

- `cp` -- polytopic Cartesian product

### Output

A list of vertices.

### Algorithm

We assume that the underlying sets are polytopic.
Then the high-dimensional set of vertices is just the Cartesian product of the
low-dimensional sets of vertices.
"""
@validate function vertices_list(cp::CartesianProduct)
    # collect low-dimensional vertices lists
    vlist1 = vertices_list(cp.X)
    vlist2 = vertices_list(cp.Y)

    # create high-dimensional vertices list
    N = eltype(cp)
    vlist = Vector{Vector{N}}()
    m = length(vlist1) * length(vlist2)
    sizehint!(vlist, m)
    for v1 in vlist1
        for v2 in vlist2
            push!(vlist, vcat(v1, v2))
        end
    end

    return vlist
end
