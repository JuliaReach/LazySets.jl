export dot_zero,
       sign_cadlag,
       remove_duplicates_sorted!,
       samedir,
       nonzero_indices,
       rectify,
       right_turn,
       is_cyclic_permutation,
       is_right_turn,
       to_negative_vector,
       _above,
       _dr,
       _up

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
    samedir(u::AbstractVector{N}, v::AbstractVector{N}) where {N<:Real}

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
                 v::AbstractVector{N}) where {N<:Real}
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
    nonzero_indices(v::AbstractVector{N}) where {N<:Real}

Return the indices in which a vector is non-zero.

### Input

- `v` -- vector

### Output

A vector of ascending indices `i` such that the vector is non-zero in dimension
`i`.
"""
function nonzero_indices(v::AbstractVector{N}) where {N<:Real}
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

function nonzero_indices(v::SparseVector{N}) where {N<:Real}
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
- `v` -- vector

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
- `Vi` -- first vector
- `Vj` -- second vector

### Output

A number indicating the direction of the difference of the given vectors.
"""
@inline function _dr(u::AbstractVector, Vi::AbstractVector, Vj::AbstractVector)
    dot(u, (Vi) - (Vj))
end

"""
    _above(u::AbstractVector, Vi::AbstractVector, Vj::AbstractVector)

Checks if the difference of the given vectors is pointing towards the given
direction.

### Input

- `u` -- direction
- `Vi` -- first vector
- `Vj` -- second vector

### Output

A boolean indicating if the difference of the given vectors is pointing
towards the given direction.
"""
@inline function _above(u::AbstractVector, Vi::AbstractVector, Vj::AbstractVector)
    _dr(u, Vi, Vj) > 0
end

"""
    to_negative_vector(v::AbstractVector{N}) where {N}

Negate a vector and convert to type `Vector`.

### Input

- `v` -- vector

### Output

A `Vector` equivalent to ``-v``.
"""
@inline function to_negative_vector(v::AbstractVector{N}) where {N}
    u = zeros(N, length(v))
    @inbounds for (i, vi) in enumerate(v)
        u[i] = -vi
    end
    return u
end

@inline function to_negative_vector(v::Vector)
    return -v
end

@inline function to_negative_vector(v::SparseVector{N}) where {N}
    u = zeros(N, length(v))
    @inbounds for (ni, i) in enumerate(v.nzind)
        u[i] = -v.nzval[ni]
    end
    return u
end

"""
    right_turn([O::AbstractVector{N}=[0, 0]], u::AbstractVector{N},
               v::AbstractVector{N}) where {N<:Real}

Compute a scalar that determines whether the acute angle defined by three 2D
points `O`, `u`, `v` in the plane is a right turn (< 180° counter-clockwise)
with respect to the center `O`.

### Input

- `O` -- (optional; default: `[0, 0]`) 2D center point
- `u` -- first 2D point
- `v` -- second 2D point

### Output

A scalar representing the rotation.
If the result is 0, the points are collinear; if it is positive, the points
constitute a positive angle of rotation around `O` from `u` to `v`; otherwise
they constitute a negative angle.

### Algorithm

The [cross product](https://en.wikipedia.org/wiki/Cross_product) is used to
determine the sense of rotation.
"""
@inline function right_turn(O::AbstractVector{N},
                            u::AbstractVector{N},
                            v::AbstractVector{N}) where {N<:Real}
    return (u[1] - O[1]) * (v[2] - O[2]) - (u[2] - O[2]) * (v[1] - O[1])
end

# version for O = origin
@inline function right_turn(u::AbstractVector{N},
                            v::AbstractVector{N}) where {N<:Real}
    return u[1] * v[2] - u[2] * v[1]
end

"""
    is_right_turn([O::AbstractVector{N}=[0, 0]], u::AbstractVector{N},
                  v::AbstractVector{N}) where {N<:Real}

Determine whether the acute angle defined by three 2D points `O`, `u`, `v`
in the plane is a right turn (< 180° counter-clockwise) with
respect to the center `O`.
Determine if the acute angle defined by two 2D vectors is a right turn (< 180°
counter-clockwise) with respect to the center `O`.

### Input

- `O` -- (optional; default: `[0, 0]`) 2D center point
- `u` -- first 2D direction
- `v` -- second 2D direction

### Output

`true` iff the vectors constitute a right turn.
"""
@inline function is_right_turn(O::AbstractVector{N},
                               u::AbstractVector{N},
                               v::AbstractVector{N}) where {N<:Real}
    return right_turn(O, u, v) >= zero(N)
end

# version for O = origin
@inline function is_right_turn(u::AbstractVector{N},
                               v::AbstractVector{N}) where {N<:Real}
    return right_turn(u, v) >= zero(N)
end
