"""
    convert(::Type{HParallelotope}, Z::AbstractZonotope)

Convert a zonotopic set of order one to a parallelotope in constraint
representation.

### Input

- `HParallelotope` -- target type
- `Z`              -- zonotopic set of order one

### Output

A parallelotope in constraint representation.

### Notes

This function requires that the list of constraints of `Z` are obtained in
the particular order returned from the `constraints_list` function of a
`Zonotope`. Hence it first converts `Z` to a `Zonotope`.
"""
function convert(::Type{HParallelotope}, Z::AbstractZonotope)
    @assert order(Z) == 1 "cannot convert a zonotope that is not of order 1 " *
                          "to a parallelotope"
    n = dim(Z)
    N = eltype(Z)

    constraints = _constraints_list_zonotope(Z)

    D = Matrix{N}(undef, n, n)
    c = Vector{N}(undef, 2n)
    j = 1
    @inbounds for i in 1:n
        D[i, :] = constraints[j].a
        c[i] = constraints[j].b
        c[i + n] = constraints[j + 1].b
        j += 2
    end
    return HParallelotope(D, c; check_consistency=false)
end
