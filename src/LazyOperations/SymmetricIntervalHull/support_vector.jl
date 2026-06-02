"""
    σ(d::AbstractVector, sih::SymmetricIntervalHull)

Return a support vector of the symmetric interval hull of a set in a given
direction.

### Input

- `d`   -- direction
- `sih` -- symmetric interval hull of a set

### Output

A support vector of the symmetric interval hull of a set in the given direction.
If the direction has norm zero, the origin is returned.

### Algorithm

For each non-zero entry in `d` we need to either look up the bound (if it has
been computed before) or compute it, in which case we store it for future
queries.
"""
@validate function σ(d::AbstractVector, sih::SymmetricIntervalHull)
    N = promote_type(eltype(d), eltype(sih))
    svec = similar(d)
    for i in eachindex(d)
        if d[i] == zero(N)
            svec[i] = zero(N)
        else
            svec[i] = sign(d[i]) * get_radius!(sih, i)
        end
    end
    return svec
end

# faster support-vector calculation for SingleEntryVector
@validate function σ(d::SingleEntryVector, sih::SymmetricIntervalHull)
    N = promote_type(eltype(d), eltype(sih))
    entry = get_radius!(sih, d.i)
    if d.v < zero(N)
        entry = -entry
    end
    return SingleEntryVector(d.i, d.n, entry)
end
