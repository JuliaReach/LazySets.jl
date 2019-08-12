export minkowski_difference, pontryagin_difference

"""
    minkowski_difference(Q::LazySet{N}, P::LazySet{N}) where {N<:Real}

Concrete Minkowski difference (geometric difference) for a pair of
lazy sets.

### Input

- `Q`         -- lazy set
- `P`         -- another lazy set that is subtracted from `Q`

### Output

An `HPolytope` that corresponds to the Minkowski difference of `Q` minus `P` if
`Q` is bounded, and an `HPolyhedron` if `Q` is unbounded.

### Notes

This function requires that the list of constraints of the lazy sets `Q` is
available and that the set `P` is bounded.
"""
function minkowski_difference(Q::LazySet{N}, P::LazySet{N}) where {N<:Real}

    @assert applicable(constraints_list, Q)  "this function " *
        "requires that the list of constraints of its first argument is applicable, but it is not; " *
        "try overapproximating with an `HPolytope` first"
    @assert isbounded(P) "this function requires that its second argument is bounded, but it is not"

    A, b = tosimplehrep(Q)
    g_QminusP = [b[i] - ρ(A[i, :], P) for i in eachindex(b)]
    if isbounded(Q)
        return HPolytope(A, g_QminusP)
    else
        return HPolyhedron(A, g_QminusP)
    end
end

const pontryagin_difference = minkowski_difference
