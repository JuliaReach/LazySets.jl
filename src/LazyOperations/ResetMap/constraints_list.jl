"""
    constraints_list(rm::ResetMap)

Return a list of constraints of a polyhedral reset map.

### Input

- `rm` -- reset map of a polyhedron

### Output

A list of constraints of the reset map.

### Notes

We assume that the underlying set `rm.X` is a polyhedron, i.e., offers a method
`constraints_list(X)`.

### Algorithm

If the set `rm.X` is hyperrectangular, we iterate through all dimensions.
For each reset we construct the corresponding (flat) constraints, and in the
other dimensions we construct the corresponding constraints of the underlying
set.

For more general sets, we fall back to `constraints_list` of a `LinearMap` of
the `A`-matrix in the affine-map view of a reset map.
Each reset dimension ``i`` is projected to zero, expressed by two constraints
for each reset dimension.
Then it remains to shift these constraints to the new value.

For instance, if the dimension ``5`` was reset to ``4``, then there will be
constraints ``x₅ ≤ 0`` and ``-x₅ ≤ 0``.
We then modify the right-hand side of these constraints to ``x₅ ≤ 4`` and
``-x₅ ≤ -4``, respectively.
"""
function constraints_list(rm::ResetMap)
    constraints = constraints_list(LinearMap(matrix(rm), set(rm)))
    N = eltype(rm)
    for (i, c) in enumerate(constraints)
        constrained_dim = find_unique_nonzero_entry(c.a)
        if constrained_dim > 0  # constrained in only one dimension
            if !haskey(rm.resets, constrained_dim)
                continue  # not a dimension we are interested in
            end
            new_value = rm.resets[constrained_dim]
            if iszero(new_value)
                @assert iszero(c.b) "expected b = 0"
                continue  # a reset to 0 needs not create a new constraint
            end
            if c.a[constrained_dim] < zero(N)
                # change sign for lower bound
                new_value = -new_value
            end
            constraints[i] = HalfSpace(c.a, new_value)
        end
    end
    return constraints
end

function constraints_list(rm::ResetMap{N,<:AbstractHyperrectangle}) where {N}
    H = rm.X
    n = dim(H)
    constraints = Vector{HalfSpace{N,SingleEntryVector{N}}}(undef, 2 * n)
    j = 1
    for i in 1:n
        ei = SingleEntryVector(i, n, one(N))
        if haskey(rm.resets, i)
            # reset dimension => add flat constraints
            v = rm.resets[i]
            constraints[j] = HalfSpace(ei, v)
            constraints[j + 1] = HalfSpace(-ei, -v)
        else
            # non-reset dimension => use the hyperrectangle's constraints
            constraints[j] = HalfSpace(ei, high(H, i))
            constraints[j + 1] = HalfSpace(-ei, -low(H, i))
        end
        j += 2
    end
    return constraints
end
