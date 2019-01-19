# default tolerance for matrix condition number (see 'isinvertible')
const DEFAULT_COND_TOL = 1e6

"""
    sign_cadlag(x::N)::N where {N<:Real}

This function works like the sign function but is ``1`` for input ``0``.

### Input

- `x` -- real scalar

### Output

``1`` if ``x ≥ 0``, ``-1`` otherwise.

### Notes

This is the sign function right-continuous at zero (see
[càdlàg function](https://en.wikipedia.org/wiki/C%C3%A0dl%C3%A0g)).
It can be used with vector-valued arguments via the dot operator.

### Examples

```jldoctest
julia> import LazySets.sign_cadlag

julia> sign_cadlag.([-0.6, 1.3, 0.0])
3-element Array{Float64,1}:
 -1.0
  1.0
  1.0
```
"""
function sign_cadlag(x::N)::N where {N<:Real}
    return x < zero(x) ? -one(x) : one(x)
end

"""
    ispermutation(u::AbstractVector{T}, v::AbstractVector{T})::Bool where T

Check that two vectors contain the same elements up to reordering.

### Input

- `u` -- first vector
- `v` -- second vector

### Output

`true` iff the vectors are identical up to reordering.

### Examples

```jldoctest
julia> LazySets.ispermutation([1, 2, 2], [2, 2, 1])
true

julia> LazySets.ispermutation([1, 2, 2], [1, 1, 2])
false

```
"""
function ispermutation(u::AbstractVector{T}, v::AbstractVector{T})::Bool where T
    if length(u) != length(v)
        return false
    end
    occurrence_map = Dict{T, Int}()
    has_duplicates = false
    for e in u
        if e ∉ v
            return false
        end
        if haskey(occurrence_map, e)
            occurrence_map[e] += 1
            has_duplicates = true
        else
            occurrence_map[e] = 1
        end
    end
    if has_duplicates
        for e in v
            if !haskey(occurrence_map, e) || occurrence_map[e] == 0
                return false
            end
            occurrence_map[e] -= 1
        end
    end
    return true
end

"""
    isinvertible(M::AbstractMatrix; [cond_tol]::Number=DEFAULT_COND_TOL)

A sufficient check of a matrix being invertible (or nonsingular).

### Input

- `M`        -- matrix
- `cond_tol` -- (optional, default: `DEFAULT_COND_TOL`) tolerance of matrix
                condition

### Output

If the result is `true`, `M` is invertible.
If the result is `false`, this function could not conclude.

### Algorithm

We check whether the
[matrix condition number](https://en.wikipedia.org/wiki/Condition_number#Matrices)
`cond(M)` is below some prescribed tolerance.
"""
function isinvertible(M::AbstractMatrix; cond_tol::Number=DEFAULT_COND_TOL)
    return cond(M) < cond_tol
end

"""
    remove_duplicates_sorted!(v::AbstractVector)

Remove duplicate entries in a sorted vector.

### Input

- `v` -- sorted vector

### Output

The input vector without duplicates.
"""
function remove_duplicates_sorted!(v::AbstractVector)
    for i in length(v)-1:-1:1
        if v[i] == v[i+1]
            splice!(v, i+1)
        end
    end
    return v
end

"""
    samedir(u::AbstractVector{N},
            v::AbstractVector{N})::Tuple{Bool, Real} where {N<:Real}

Check whether two vectors point in the same direction.

### Input

- `u` -- first vector
- `v` -- second vector

### Output

`(true, k)` iff the vectors are identical up to a positive scaling factor `k`,
and `(false, 0)` otherwise.


### Examples

```jldoctest
julia> LazySets.samedir([1, 2, 3], [2, 4, 6])
(true, 0.5)

julia> LazySets.samedir([1, 2, 3], [3, 2, 1])
(false, 0)

julia> LazySets.samedir([1, 2, 3], [-1, -2, -3])
(false, 0)

```
"""
function samedir(u::AbstractVector{N},
                 v::AbstractVector{N}
                )::Tuple{Bool, Real} where {N<:Real}
    @assert length(u) == length(v) "wrong dimension"
    no_factor = true
    factor = 0
    @inbounds for i in 1:length(u)
        if u[i] == 0
            if v[i] != 0
                return (false, 0)
            end
            continue
        elseif v[i] == 0
            return (false, 0)
        end
        if no_factor
            no_factor = false
            factor = u[i] / v[i]
            if factor < 0
                return (false, 0)
            end
        elseif factor != u[i] / v[i]
            return (false, 0)
        end
    end
    if no_factor
        # both vectors are zero
        return (true, 0)
    end
    return (true, factor)
end

"""
    nonzero_indices(v::AbstractVector{N})::Vector{Int} where {N<:Real}

Return the indices in which a vector is non-zero.

### Input

- `v` -- vector

### Output

A vector of ascending indices `i` such that the vector is non-zero in dimension
`i`.
"""
function nonzero_indices(v::AbstractVector{N})::Vector{Int} where {N<:Real}
    n = length(v)
    res = Vector{Int}()
    sizehint!(res, n)
    for i in 1:n
        if v[i] != zero(N)
            push!(res, i)
        end
    end
    return res
end

function nonzero_indices(v::SparseVector{N})::Vector{Int} where {N<:Real}
    return x.nzind
end

"""
    reseed(rng::AbstractRNG, seed::Union{Int, Nothing})::AbstractRNG

Reset the RNG seed if the seed argument is a number.

### Input

- `rng`  -- random number generator
- `seed` -- seed for reseeding

### Output

The input RNG if the seed is `nothing`, and a reseeded RNG otherwise.
"""
function reseed(rng::AbstractRNG, seed::Union{Int, Nothing})::AbstractRNG
    if seed != nothing
        @static if VERSION < v"0.7-"
            return Random.srand(rng, seed)
        else
            return Random.seed!(rng, seed)
        end
    end
    return rng
end

"""
    cross_product(M::AbstractMatrix{N}) where {N<:Real}

Compute the high-dimensional cross product of ``n-1`` ``n``-dimensional vectors.

### Input

- `M` -- ``n × n - 1``-dimensional matrix

### Output

A vector.

### Algorithm

The cross product is defined as follows:

```math
\\left[ \\dots, (-1)^{n+1} \\det(M^{[i]}), \\dots \\right]^T
```
where ``M^{[i]}`` is defined as ``M`` with the ``i``-th row removed.
See *Althoff, Stursberg, Buss: Computing Reachable Sets of Hybrid Systems Using
a Combination of Zonotopes and Polytopes. 2009.*
"""
function cross_product(M::AbstractMatrix{N})::Vector{N} where {N<:Real}
    n = size(M, 1)
    @assert 1 < n == size(M, 2) + 1 "the matrix must be n x (n-1) dimensional"

    v = Vector{N}(undef, n)
    minus = false
    for i in 1:n
        Mi = view(M, 1:n .!= i, :)  # remove i-th row
        d = det(Mi)
        if minus
            v[i] = -d
            minus = false
        else
            v[i] = d
            minus = true
        end
    end
    return v
end

"""
    StrictlyIncreasingIndices

Iterator over the vectors of `m` strictly increasing indices from 1 to `n`.

### Fields

- `n` -- size of the index domain
- `m` -- number of indices to choose (resp. length of the vectors)

### Notes

The vectors are modified in-place.

The iterator ranges over ``\\binom{n}{m}`` (`n` choose `m`) possible vectors.

This implementation results in a lexicographic order with the last index growing
first.

### Examples

```jldoctest
julia> for v in LazySets.StrictlyIncreasingIndices(4, 2)
           println(v)
       end
[1, 2]
[1, 3]
[1, 4]
[2, 3]
[2, 4]
[3, 4]
```
"""
struct StrictlyIncreasingIndices
    n::Int
    m::Int

    function StrictlyIncreasingIndices(n::Int, m::Int)
        @assert n >= m > 0 "require n >= m > 0"
        new(n, m)
    end
end

Base.eltype(::Type{StrictlyIncreasingIndices}) = Vector{Int}
Base.length(sii::StrictlyIncreasingIndices) = binomial(sii.n, sii.m)

@static if VERSION < v"0.7-"
@eval begin

# returns [1, 2, ..., m-2, m-1, m-1]
function Base.start(sii::StrictlyIncreasingIndices)
    v = [1:sii.m;]
    if sii.n > sii.m
        v[end] -= 1
    end
    return v
end

function Base.next(sii::StrictlyIncreasingIndices, state)
    v = state
    i = sii.m
    diff = sii.n
    if i == diff
        return (v, nothing)
    end
    while v[i] == diff
        i -= 1
        diff -= 1
    end
    # update vector
    v[i] += 1
    for j in i+1:sii.m
        v[j] = v[j-1] + 1
    end
    # detect termination: first index has maximum value
    if i == 1 && v[1] == (sii.n - sii.m + 1)
        return (v, nothing)
    end
    return (v, v)
end

Base.done(sii::StrictlyIncreasingIndices, state) = state == nothing

end # @eval
else
@eval begin

# initialization
function Base.iterate(sii::StrictlyIncreasingIndices)
    v = [1:sii.m;]
    return (v, v)
end

# normal iteration
function Base.iterate(sii::StrictlyIncreasingIndices, state::AbstractVector{Int})
    v = state
    i = sii.m
    diff = sii.n
    if i == diff
        return nothing
    end
    while v[i] == diff
        i -= 1
        diff -= 1
    end
    # update vector
    v[i] += 1
    for j in i+1:sii.m
        v[j] = v[j-1] + 1
    end
    # detect termination: first index has maximum value
    if i == 1 && v[1] == (sii.n - sii.m + 1)
        return (v, nothing)
    end
    return (v, v)
end

# termination
function Base.iterate(sii::StrictlyIncreasingIndices, state::Nothing)
    return nothing
end

end # @eval
end # if
