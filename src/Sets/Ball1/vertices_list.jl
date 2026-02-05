"""
# Extended help

    vertices_list(B::Ball1)

### Notes

In ``n`` dimensions there are ``2n`` vertices (unless the radius is 0).
"""
@validate function vertices_list(B::Ball1)
    # fast evaluation if B has radius 0
    if iszero(B.radius)
        return [B.center]
    end
    VN = _vector_type(B)
    vertices = Vector{VN}(undef, 2 * dim(B))
    j = 0
    v = copy(B.center)
    @inbounds for i in 1:dim(B)
        ci = v[i]
        v[i] += B.radius
        j += 1
        vertices[j] = copy(v)
        v[i] = ci - B.radius
        j += 1
        vertices[j] = copy(v)
        v[i] = ci  # restore old value
    end
    return vertices
end

_vector_type(B::Ball1{N,VN}) where {N,VN} = VN
