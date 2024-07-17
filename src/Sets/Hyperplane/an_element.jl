"""
    an_element(H::Hyperplane)

Return some element of a hyperplane.

### Input

- `H` -- hyperplane

### Output

An element on the hyperplane.
"""
function an_element(H::Hyperplane)
    return _an_element_helper_hyperplane(H.a, H.b)
end

"""
    _an_element_helper_hyperplane(a::AbstractVector{N}, b,
                                  [nonzero_entry_a]::Int) where {N}

Helper function that computes an element on a hyperplane ``a⋅x = b``.

### Input

- `a`               -- normal direction
- `b`               -- constraint
- `nonzero_entry_a` -- (optional, default: computes the first index) index `i`
                       such that `a[i]` is different from 0

### Output

An element on a hyperplane.

### Algorithm

We compute the point on the hyperplane as follows:
- We already found a nonzero entry of ``a`` in dimension, say, ``i``.
- We set ``x[i] = b / a[i]``.
- We set ``x[j] = 0`` for all ``j ≠ i``.
"""
@inline function _an_element_helper_hyperplane(a::AbstractVector{N}, b,
                                               nonzero_entry_a::Int=findfirst(!iszero, a)) where {N}
    x = zeros(N, length(a))
    x[nonzero_entry_a] = b / a[nonzero_entry_a]
    return x
end
