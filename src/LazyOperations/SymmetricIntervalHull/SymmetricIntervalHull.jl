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

The convenience alias `⊡` can be typed by `\\boxdot<tab>`.
"""
struct SymmetricIntervalHull{N,S<:LazySet{N}} <: AbstractHyperrectangle{N}
    X::S
    cache::Vector{N}

    # default constructor that initializes cache
    function SymmetricIntervalHull(X::S;
                                   check_boundedness::Bool=true) where {N,S<:LazySet{N}}
        @assert !check_boundedness || isbounded(X) "the symmetric interval hull is only defined " *
                                                   "*for bounded sets"

        # fill cache with default value -1 (actual bounds cannot be negative)
        cache = fill(-one(N), dim(X))

        return new{N,S}(X, cache)
    end
end

const ⊡ = SymmetricIntervalHull

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
@validate function radius_hyperrectangle(sih::SymmetricIntervalHull, i::Int)
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

include("center.jl")
include("dim.jl")
include("isoperationtype.jl")
include("support_function.jl")
include("support_vector.jl")
