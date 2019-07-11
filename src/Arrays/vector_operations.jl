export dot_zero,
       sign_cadlag,
       ispermutation,
       remove_duplicates_sorted!,
       samedir,
       nonzero_indices,
       rectify,
       is_cyclic_permutation

"""
    dot_zero(x::AbstractVector{N}, y::AbstractVector{N}) where{N<:Real}

Dot product with preference for zero value in the presence of infinity values.

### Input

- `x` -- first vector
- `y` -- second vector

### Output

The dot product of `x` and `y`, but with the rule that `0 * Inf == 0`.
"""
function dot_zero(x::AbstractVector{N}, y::AbstractVector{N}) where{N<:Real}
    res = zero(N)
    for i in 1:length(x)
        if !iszero(x[i]) && !iszero(y[i])
            res += x[i] * y[i]
        end
    end
    return res
end

"""
    ispermutation(u::AbstractVector{T}, v::AbstractVector)::Bool where {T}

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
function ispermutation(u::AbstractVector{T}, v::AbstractVector)::Bool where {T}
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
julia> using LazySets: samedir

julia> samedir([1, 2, 3], [2, 4, 6])
(true, 0.5)

julia> samedir([1, 2, 3], [3, 2, 1])
(false, 0)

julia> samedir([1, 2, 3], [-1, -2, -3])
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
    return v.nzind
end

"""
    rectify(x::AbstractVector{N}) where {N<:Real}

Rectify a vector, i.e., take the element-wise maximum with zero.

### Input

- `x` -- vector

### Output

A copy of the vector where each negative entry is replaced by zero.
"""
function rectify(x::AbstractVector{N}) where {N<:Real}
    return map(xi -> max(xi, zero(N)), x)
end

"""
    is_cyclic_permutation(candidate::AbstractVector, paragon::AbstractVector)

Checks if the elements in `candidate` are a cyclic permutation of the elements
in `paragon`.

### Input

- `candidate` -- candidate vector
- `paragon`   -- paragon vector

### Output

A boolean indicating if the elements of `candidate` are in the same order as in
`paragon` or any of its cyclic permutations.
"""
function is_cyclic_permutation(candidate::AbstractVector,
                               paragon::AbstractVector)
    m = length(candidate)
    if length(paragon) != m
        return false
    end
    return any(candidate == circshift(paragon, i) for i in 0:m-1)
end

"""
    _up(u::AbstractVector, v::AbstractVector)

Checks if the given vector is pointing towards the given direction.

### Input

- `u` -- direction
- `v`   -- vector

### Output

A boolean indicating if the vector is pointing towards the direction.
"""
@inline function _up(u::AbstractVector, v::AbstractVector)
    dot(u, v) > 0
end

"""
    _dr(u::AbstractVector, Vi::AbstractVector, Vj::AbstractVector)

Returns the direction of the difference of the given vectors.

### Input

- `u` -- direction
- `Vi`   -- first vector
- `Vj` -- second vector

### Output

A float indicating the direction of the difference of the given vectors.
"""
@inline function _dr(u::AbstractVector, Vi::AbstractVector, Vj::AbstractVector)
    (dot(u, (Vi) - (Vj)))
end

"""
    _above(u::AbstractVector, Vi::AbstractVector, Vj::AbstractVector)

Checks if the difference of the given vectors are pointing towards the given
direction.

### Input

- `u` -- direction
- `Vi`   -- first vector
- `Vj` -- second vector

### Output

A boolean indicating if the difference of the given vectors are pointing
towards the given direction.
"""
@inline function _above(u::AbstractVector, Vi::AbstractVector, Vj::AbstractVector)
    (_dr(u, Vi, Vj) > 0)
end
