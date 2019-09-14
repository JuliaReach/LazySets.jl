export reflect

"""
    reflect(P::LazySet{N}) where {N<:Real}

Concrete reflection of a convex set `P`, resulting in the reflected set `-P`.

### Note

This function requires that the list of constraints of the set `P` is
available.

This function can be used to implement the alternative definition of the
Minkowski Difference, which writes as
```math
A ⊖ B = \\{a − b | a ∈ A, b ∈ B} = A ⊕ (-B)
P = \\{z ∈ ℝⁿ: sᵢᵀz ≤ rᵢ, i=1,...,N\\}.
```
by calling `minkowski_sum(A, reflect(B))`.
"""
function reflect(P::LazySet{N}) where {N<:Real}

    @assert applicable(constraints_list, P)  "this function " *
        "requires that the list of constraints of its first argument is applicable, but it is not; " *
        "if it is bounded, try overapproximating with an `HPolytope` first"

    F,g = tosimplehrep(P)
    if isbounded(P)
        return HPolytope(-F, g)
    else
        return HPolyhedron(-F, g)
    end
end
