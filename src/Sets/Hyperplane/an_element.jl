"""
# Extended help

    an_element(H::Hyperplane)

### Algorithm

We compute a point on the hyperplane ``a⋅x = b`` as follows:
- We first find a nonzero entry of ``a`` in dimension, say, ``i``.
- We set ``x[i] = b / a[i]``.
- We set ``x[j] = 0`` for all ``j ≠ i``.
"""
function an_element(H::Hyperplane)
    return _an_element_helper_hyperplane(H.a, H.b)
end

@inline function _an_element_helper_hyperplane(a::AbstractVector{N}, b,
                                               nonzero_entry_a::Int=findfirst(!iszero, a)) where {N}
    x = zeros(N, length(a))
    x[nonzero_entry_a] = b / a[nonzero_entry_a]
    return x
end
