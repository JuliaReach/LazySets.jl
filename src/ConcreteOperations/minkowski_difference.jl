export minkowski_difference, pontryagin_difference

"""
    minkowski_difference(P::LazySet, Q::LazySet)

Concrete Minkowski difference (geometric difference) of a polytopic set and a
compact convex set.

### Input

- `P` -- polytopic set
- `Q` -- compact convex set that is subtracted from `P`

### Output

An `HPolytope` that corresponds to the Minkowski difference of `P` minus `Q` if
`P` is bounded, and an `HPolyhedron` if `P` is unbounded.

### Notes

This implementation requires that the set `P` is polyhedral and that the set `Q`
is bounded.

### Algorithm

This method implements Theorem 2.3 in [1]:

Suppose ``P`` is a polyhedron
```math
P = \\{z ∈ ℝ^n: sᵢᵀz ≤ rᵢ,~i = 1, …, N\\}.
```
where ``sᵢ ∈ ℝ^n, sᵢ ≠ 0``, and ``rᵢ ∈ ℝ``.
Assume ``ρ(sᵢ,Q)`` is defined for ``i = 1, …, N``.
Then the Minkowski difference is

```math
\\{z ∈ ℝ^n: sᵢᵀz ≤ rᵢ - ρ(sᵢ,Q),~i = 1, …, N\\}.
```

[1] Ilya Kolmanovsky and Elmer G. Gilbert (1997). *Theory and computation
of disturbance invariant sets for discrete-time linear systems.*
[Mathematical Problems in Engineering Volume 4, Issue 4, Pages
317-367.](http://dx.doi.org/10.1155/S1024123X98000866)
"""
function minkowski_difference(P::LazySet, Q::LazySet)

    @assert is_polyhedral(P)  "this implementation requires that the first" *
        "argument is polyhedral; try overapproximating with an `HPolyhedron`"
    @assert isbounded(Q) "this implementation requires that the second " *
        "argument is bounded, but it is not"

    A, b = tosimplehrep(P)
    g_PminusQ = [b[i] - ρ(A[i, :], Q) for i in eachindex(b)]
    if isbounded(P)
        return HPolytope(A, g_PminusQ)
    else
        return HPolyhedron(A, g_PminusQ)
    end
end

"""
    pontryagin_difference(X::LazySet, Y::LazySet)

An alias for the function `minkowski_difference`.

### Notes

There is some inconsistency in the literature regarding the naming conventions.
In this library, both the names *Minkowski difference* and
*Pontryagin difference* refer to the geometric difference of two sets.
Mathematically:

```math
    X ⊖ Y = \\{z ∈ ℝ^n: z + v ∈ X  ~∀~v ∈ Y\\}
```
"""
const pontryagin_difference = minkowski_difference

# concrete Minkowski difference with singleton
minkowski_difference(X::LazySet, S::AbstractSingleton) = translate(X, -element(S))
minkowski_difference(X::LazySet, ::ZeroSet) = X
