export SymmetricIntervalHull, ⊡

"""
    SymmetricIntervalHull{N, S<:LazySet{N}} <: AbstractHyperrectangle{N}

Type that represents the symmetric interval hull of a compact set.

### Fields

- `X`     -- compact set
- `cache` -- partial storage of already computed bounds, organized as mapping
             from the dimension to the bound value

### Notes

The symmetric interval hull can be computed with ``2n`` support-function queries
(of unit vectors), where ``n`` is the dimension of the wrapped set (i.e., two
queries per dimension).
When asking for the support vector (or support function) in a direction ``d``,
one needs ``2k`` such queries, where ``k`` is the number of non-zero entries in
``d``.

However, if one asks for many support vectors (or support-function evaluations)
in a loop, the number of computations may exceed ``2n``.
To be most efficient in such cases, this type stores the intermediately computed
bounds in the `cache` field.

The set `X` must be bounded. The flag `check_boundedness` (which defaults to
`true`) can be used to elide the boundedness check in the inner constructor.
Misuse of this flag can result in incorrect behavior.

The symmetric interval hull of a set is a hyperrectangle centered in the origin,
which in particular is convex.

An alias for this function is `⊡`.
"""
struct SymmetricIntervalHull{N,S<:LazySet{N}} <: AbstractHyperrectangle{N}
    X::S
    cache::Vector{N}

    # default constructor that initializes cache
    function SymmetricIntervalHull(X::S;
                                   check_boundedness::Bool=true) where {N,S<:LazySet{N}}
        @assert !check_boundedness || isboundedtype(typeof(X)) ||
                isbounded(X) "the symmetric interval hull is only defined " *
                             "for bounded sets"

        # fill cache with default value -1 (actual bounds cannot be negative)
        cache = fill(-one(N), dim(X))

        return new{N,S}(X, cache)
    end
end

"""
    ⊡

Alias for `SymmetricIntervalHull`.
"""
const ⊡ = SymmetricIntervalHull

isoperationtype(::Type{<:SymmetricIntervalHull}) = true

"""
    SymmetricIntervalHull(∅::EmptySet)

The symmetric interval hull of an empty set.

### Input

- `∅` -- an empty set

### Output

The empty set because it is absorbing for the symmetric interval hull.
"""
SymmetricIntervalHull(∅::EmptySet) = ∅

"""
    radius_hyperrectangle(sih::SymmetricIntervalHull, i::Int)

Return the box radius of the symmetric interval hull of a set in a given
dimension.

### Input

- `sih` -- symmetric interval hull of a set
- `i`   -- dimension of interest

### Output

The radius in the given dimension.

### Notes

If the radius was computed before, this is just a look-up. Otherwise it is
computed.
"""
function radius_hyperrectangle(sih::SymmetricIntervalHull, i::Int)
    return get_radius!(sih, i)
end

"""
    radius_hyperrectangle(sih::SymmetricIntervalHull)

Return the box radius of the symmetric interval hull of a set in every
dimension.

### Input

- `sih` -- symmetric interval hull of a set

### Output

The box radius of the symmetric interval hull of a set.

### Notes

This function computes the symmetric interval hull explicitly.
"""
function radius_hyperrectangle(sih::SymmetricIntervalHull)
    for i in 1:dim(sih)
        get_radius!(sih, i)
    end
    return sih.cache
end

"""
    center(sih::SymmetricIntervalHull, i::Int)

Return the center along a given dimension of the symmetric interval hull of a
set.

### Input

- `sih` -- symmetric interval hull of a set
- `i`   -- dimension of interest

### Output

The center along a given dimension of the symmetric interval hull of a set.
"""
@inline function center(sih::SymmetricIntervalHull, i::Int)
    @boundscheck _check_bounds(sih, i)
    N = eltype(sih)
    return zero(N)
end

"""
    center(sih::SymmetricIntervalHull)

Return the center of the symmetric interval hull of a set.

### Input

- `sih` -- symmetric interval hull of a set

### Output

The origin.
"""
function center(sih::SymmetricIntervalHull)
    N = eltype(sih)
    return zeros(N, dim(sih))
end

"""
    dim(sih::SymmetricIntervalHull)

Return the dimension of the symmetric interval hull of a set.

### Input

- `sih` -- symmetric interval hull of a set

### Output

The ambient dimension of the symmetric interval hull of a set.
"""
function dim(sih::SymmetricIntervalHull)
    return dim(sih.X)
end

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

# faster support-function evaluation for SingleEntryVector
@validate function ρ(d::SingleEntryVector, sih::SymmetricIntervalHull)
    return abs(d.v) * get_radius!(sih, d.i)
end

"""
    get_radius!(sih::SymmetricIntervalHull, i::Int)

Compute the radius of the symmetric interval hull of a set in a given dimension.

### Input

- `sih` -- symmetric interval hull of a set
- `i`   -- dimension in which the radius should be computed

### Output

The radius of the symmetric interval hull of a set in a given dimension.

### Algorithm

We ask for the `extrema` of the underlying set in dimension `i`.
"""
function get_radius!(sih::SymmetricIntervalHull, i::Int)
    N = eltype(sih)
    if sih.cache[i] == -one(N)
        # compute the radius
        l, h = extrema(sih.X, i)
        r = max(h, abs(l))
        sih.cache[i] = r
        return r
    end
    return sih.cache[i]
end
