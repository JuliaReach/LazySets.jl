export SymmetricIntervalHull,
       an_element

"""
    SymmetricIntervalHull{N<:Real, S<:LazySet{N}} <: LazySet{N}

Type that represents the symmetric interval hull of a convex set.

### Fields

- `X`      -- convex set
- `radius` -- partial storage of already computed bounds
- `known`  -- bit array that marks already computed bounds

### Notes

The symmetric interval hull can be computed with ``2n`` support vector queries
of unit vectors, where ``n`` is the dimension of the wrapped set (i.e., two
queries per dimension).
When asking for the support vector for a direction ``d``, one needs ``2k`` such
queries, where ``k`` is the number of non-zero entries in ``d``.

However, if one asks for many support vectors in a loop, the number of
computations may exceed ``2n``.
To be most efficient in such cases, this type stores the intermediately computed
bounds in the `radius` field and maintains a `known` bit array to mark already
computed bounds.
"""
struct SymmetricIntervalHull{N<:Real, S<:LazySet{N}} <: LazySet{N}
    X::S
    radius::Vector{N}
    known::Vector{Bool}

    function SymmetricIntervalHull{N, S}(X::S) where {S<:LazySet{N}} where {N<:Real}
        n = dim(X)
        return new{N, S}(X, zeros(N, n), falses(n))
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
        elseif !sih.known[i]
            # compute support vector in dimension i
            right_bound = σ(sparsevec([i], [one_N], len), sih.X)
            left_bound = σ(sparsevec([i], [-one_N], len), sih.X)
            sih.radius[i] = max(right_bound[i], abs(left_bound[i]))
            sih.known[i] = true
        end
        svec[i] = sign(d[i]) * sih.radius[i]
    end
    return svec
end

"""
    an_element(sih::SymmetricIntervalHull{N, S}
              )::Vector{N} where {N<:Real, S<:LazySet{N}}

Return some element of a symmetric interval hull.

### Input

- `sih` -- symmetric interval hull

### Output

The origin.
"""
function an_element(sih::SymmetricIntervalHull{N, S}
                   )::Vector{N} where {N<:Real, S<:LazySet{N}}
    return zeros(N, dim(sih))
end
