export minkowski_difference, pontryagin_difference

"""
    minkowski_difference(P::LazySet{N}, Q::LazySet{N}) where {N<:Real}

Concrete Minkowski difference (geometric difference) for a pair of
convex sets.

### Input

- `P` -- polytopic set
- `Q` -- compact convex set that is subtracted from `P`

### Output

An `HPolytope` that corresponds to the Minkowski difference of `P` minus `Q` if
`P` is bounded, and an `HPolyhedron` if `P` is unbounded.

### Notes

This function requires that the list of constraints of the lazy sets `P` is
available and that the set `Q` is bounded.

### Algorithm

This function implements [Theorem 2.2], that we state next.
Suppose ``U`` is compact and convex and ``U ~ V ≠ Ø``.
Then
i) ``hᵥ(η)`` and ``hᵤ(η)`` are defined for all ``η ∈ ℝⁿ``;
ii)
```math
U ~ V = {z ∈ ℝⁿ: ηᵀz ≤ hᵤ(η) - hᵥ(η) ∀η ∈ ℝⁿ }.
```
where ``~`` is defined as ``U ~ V = {z ∈ ℝⁿ: z + v ∈ U  ∀ v ∈ V}`` and is called
the *Minkowski difference* (also referenced as *Pontryagin difference*)
and often denoted by the operator ⊖.

[1] *Ilya Kolmanovsky and Elmer G. Gilbert (1997). Theory and computation
of disturbance invariant sets for discrete-time linear systems.
Mathematical Problems in Engineering Volume 4, Issue 4, Pages 317-367.
http://dx.doi.org/10.1155/S1024123X98000866*
"""
function minkowski_difference(P::LazySet{N}, Q::LazySet{N}) where {N<:Real}

    @assert applicable(constraints_list, P)  "this function " *
        "requires that the list of constraints of its first argument is applicable, but it is not; " *
        "if it is bounded, try overapproximating with an `HPolytope` first"
    @assert isbounded(Q) "this function requires that its second argument is bounded, but it is not"

    A, b = tosimplehrep(P)
    g_PminusQ = [b[i] - ρ(A[i, :], Q) for i in eachindex(b)]
    if isbounded(P)
        return HPolytope(A, g_PminusQ)
    else
        return HPolyhedron(A, g_PminusQ)
    end
end


"""
    pontryagin_difference(P::LazySet{N}, Q::LazySet{N}) where {N<:Real}

An alias for the function `minkowski_difference`.

### Notes

Due to inconsistent naming conventions, both the name Minkowski Difference and
Pontryagin Difference are used to refer to the geometric difference of two sets.

"""
const pontryagin_difference = minkowski_difference
