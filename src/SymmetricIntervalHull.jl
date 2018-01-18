export SymmetricIntervalHull,
       an_element

"""
    SymmetricIntervalHull{N<:Real, S<:LazySet{N}} <: AbstractHyperrectangle{N}

Type that represents the symmetric interval hull of a convex set.

### Fields

- `X`     -- convex set
- `cache` -- partial storage of already computed bounds, organized as mapping
    from dimension to tuples `(bound, valid)`, where `valid` is a flag
    indicating if the `bound` entry has been computed

### Notes

The symmetric interval hull can be computed with ``2n`` support vector queries
of unit vectors, where ``n`` is the dimension of the wrapped set (i.e., two
queries per dimension).
When asking for the support vector for a direction ``d``, one needs ``2k`` such
queries, where ``k`` is the number of non-zero entries in ``d``.

However, if one asks for many support vectors in a loop, the number of
computations may exceed ``2n``.
To be most efficient in such cases, this type stores the intermediately computed
bounds in the `cache` field.
"""
struct SymmetricIntervalHull{N<:Real, S<:LazySet{N}} <: AbstractHyperrectangle{N}
    X::S
    cache::Vector{Tuple{N, Bool}}

    function SymmetricIntervalHull{N, S}(X::S) where {S<:LazySet{N}} where {N<:Real}
        cache = Vector{Tuple{N, Bool}}(dim(X))
        zero_N = zero(N)
        for i in 1:length(cache)
            cache[i] = (zero_N, false)
        end
        return new{N, S}(X, cache)
    end
end
# type-less convenience constructor
SymmetricIntervalHull(X::S) where {S<:LazySet{N}} where {N<:Real} =
    SymmetricIntervalHull{N, S}(X)

"""
    SymmetricIntervalHull(∅::EmptySet)::EmptySet

The symmetric interval hull of an empty set.

### Input

- `∅` -- an empty set

### Output

The empty set because it is absorbing for the symmetric interval hull.
"""
SymmetricIntervalHull(∅::EmptySet)::EmptySet = ∅


# --- AbstractHyperrectangle interface functions ---


"""
    radius_hyperrectangle(H::SymmetricIntervalHull{N}, i::Int)::N where {N<:Real}

Return the box radius of a symmetric interval hull of a convex set in a given
dimension.

### Input

- `H` -- symmetric interval hull of a convex set

### Output

The radius in the given dimension.
If it was computed before, this is just a look-up, otherwise it requires two
support vector computations.
"""
function radius_hyperrectangle(H::SymmetricIntervalHull{N},
                               i::Int)::N where {N<:Real}
    if !sih.cache[i][2]
        compute_radius(sih, i)
    end
    return sih.cache[i][1]
end

"""
    radius_hyperrectangle(sih::SymmetricIntervalHull{N})::Vector{N} where {N<:Real}

Return the box radius of a symmetric interval hull of a convex set in every
dimension.

### Input

- `sih` -- symmetric interval hull of a convex set

### Output

The box radius of the symmetric interval hull of a convex set.

### Notes

This function computes the symmetric interval hull explicitly.
"""
function radius_hyperrectangle(sih::SymmetricIntervalHull{N}
                              )::Vector{N} where {N<:Real}
    dim = dim(sih)
    one_N = one(N)
    for i in 1:dim
        if !sih.cache[i][2]
            compute_radius(sih, i, dim, one_N)
        end
    end
    return sih.cache[i][1]
end


# --- AbstractPointSymmetric interface functions ---


"""
    center(sih::SymmetricIntervalHull{N})::Vector{N} where {N<:Real}

Return the center of a symmetric interval hull of a convex set.

### Input

- `sih` -- symmetric interval hull of a convex set

### Output

The origin.
"""
function center(sih::SymmetricIntervalHull{N})::Vector{N} where {N<:Real}
    return zeros(N, dim(sih))
end


# --- LazySet interface functions ---


"""
    dim(sih::SymmetricIntervalHull)::Int

Return the dimension of a symmetric interval hull of a convex set.

### Input

- `sih` -- symmetric interval hull of a convex set

### Output

The ambient dimension of the symmetric interval hull of a convex set.
"""
function dim(sih::SymmetricIntervalHull)::Int
    return dim(sih.X)
end

"""
    σ(d::AbstractVector{N}, sih::SymmetricIntervalHull
     )::AbstractVector{<:Real} where {N<:Real}

Return the support vector of a symmetric interval hull of a convex set in a
given direction.

### Input

- `d`   -- direction
- `sih` -- symmetric interval hull of a convex set

### Output

The support vector of the symmetric interval hull of a convex set in the given
direction.
If the direction has norm zero, the origin is returned.

### Algorithm

For each non-zero entry in `d` we need to either look up the bound (if it has
been computed before) or compute it, in which case we store it for future
queries.
One such computation just asks for the support vector of the underlying set for
both the positive and negative unit vector in the respective dimension.
"""
function σ(d::AbstractVector{N}, sih::SymmetricIntervalHull
          )::AbstractVector{<:Real} where {N<:Real}
    len = length(d)
    @assert len == dim(sih) "cannot compute the support vector of a " *
        "$(dim(sih))-dimensional set along a vector of length $(length(d))"

    svec = similar(d)
    zero_N = zero(N)
    one_N = one(N)
    for i in eachindex(d)
        if d[i] == zero_N
            svec[i] = zero_N
        else
            if !sih.cache[i][2]
                compute_radius(sih, i, len, one_N)
            end
            svec[i] = sign(d[i]) * sih.cache[i][1]
        end
    end
    return svec
end


# --- SymmetricIntervalHull functions ---


"""
    compute_radius(sih::SymmetricIntervalHull{N, <:LazySet{N}},
                   i::Int,
                   dimension::Int=dim(sih),
                   one_N::N=one(N))::Void where {N<:Real}

Compute the radius of a symmetric interval hull of a convex in a given
dimension.

### Input

- `sih`       -- symmetric interval hull of a convex set
- `i`         -- dimension
- `dimension` -- (optional, default: `dim(sih)`) set dimension
- `one_N`     -- (optional, default: `one(N)`) numeric representation of `1`

### Output

Nothing.

### Algorithm

We ask for the support vector of the underlying set for both the positive and
negative unit vector in the dimension `i`.
"""
function compute_radius(sih::SymmetricIntervalHull{N, <:LazySet{N}},
                        i::Int,
                        dimension::Int=dim(sih),
                        one_N::N=one(N))::Void where {N<:Real}
    right_bound = σ(sparsevec([i], [one_N], dimension), sih.X)
    left_bound = σ(sparsevec([i], [-one_N], dimension), sih.X)
    sih.cache[i] = (max(right_bound[i], abs(left_bound[i])), true)
    return nothing
end
