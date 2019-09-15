export reflect

"""
    reflect(P::LazySet{N}) where {N<:Real}

Concrete reflection of a convex set `P`, resulting in the reflected set `-P`.

### Note

This function requires that the list of constraints of the set `P` is
available, i.e. such that it can be written as
``P = \\{z ∈ ℝⁿ: ⋂ sᵢᵀz ≤ rᵢ, i = 1, ..., N\\}.``

This function can be used to implement the alternative definition of the
Minkowski Difference, which writes as
```math
A ⊖ B = \\{a − b | a ∈ A, b ∈ B\\} = A ⊕ (-B)
```
by calling `minkowski_sum(A, reflect(B))`.
"""
function reflect(P::LazySet)

    @assert applicable(constraints_list, P)  "this function " *
        "requires that the list of constraints is available, but it is not; " *
        "if the set is bounded, try overapproximating with an `HPolytope` first"

    F,g = tosimplehrep(P)
    T = isbounded(P) ? HPolytope : HPolyhedron
    return T(-F, g)
end
