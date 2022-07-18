export dot_zero,
       remove_duplicates_sorted!,
       samedir,
       ismultiple,
       nonzero_indices,
       right_turn,
       is_cyclic_permutation,
       is_right_turn,
       to_negative_vector,
       _above,
       _dr,
       _up,
       distance,
       append_zeros,
       prepend_zeros

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
    samedir(u::AbstractVector{<:Real}, v::AbstractVector{<:Real})

Check whether two vectors point in the same direction.

### Input

- `u` -- first vector
- `v` -- second vector

### Output

`(true, k)` iff the vectors are identical up to a positive scaling factor `k`
such that `u = k * v`, and `(false, 0)` otherwise.


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
function samedir(u::AbstractVector{<:Real}, v::AbstractVector{<:Real})
    return _ismultiple(u, v; allow_negative=false)
end

"""
    ismultiple(u::AbstractVector{<:Real}, v::AbstractVector{<:Real})

Check whether two vectors are linearly dependent.

### Input

- `u` -- first vector
- `v` -- second vector

### Output

`(true, k)` iff the vectors are identical up to a scaling factor `k ≠ 0` such
that `u = k * v`, and `(false, 0)` otherwise.


### Examples

```jldoctest
julia> using LazySets: ismultiple

julia> ismultiple([1, 2, 3], [2, 4, 6])
(true, 0.5)

julia> ismultiple([1, 2, 3], [3, 2, 1])
(false, 0)

julia> ismultiple([1, 2, 3], [-1, -2, -3])
(true, -1.0)

```
"""
function ismultiple(u::AbstractVector{<:Real}, v::AbstractVector{<:Real})
    return _ismultiple(u, v; allow_negative=true)
end

function _ismultiple(u::AbstractVector, v::AbstractVector; allow_negative::Bool)
    @assert length(u) == length(v) "wrong dimension"
    no_factor = true
    factor = 0
    @inbounds for i in 1:length(u)
        if isapproxzero(u[i])
            if !isapproxzero(v[i])
                return (false, 0)
            end
            continue
        elseif isapproxzero(v[i])
            return (false, 0)
        end
        if no_factor
            no_factor = false
            factor = u[i] / v[i]
            if !allow_negative && factor < 0
                return (false, 0)
            end
        elseif !_isapprox(factor, u[i] / v[i])
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
    return _geq(right_turn(O, u, v), zero(N))
end

# version for O = origin
@inline function is_right_turn(u::AbstractVector{N},
                               v::AbstractVector{N}) where {N<:Real}
    return _geq(right_turn(u, v), zero(N))
end

"""
    distance(x::AbstractVector, y::AbstractVector; [p]::Real=2.0)

Compute the distance between two vectors with respect to the given `p`-norm,
computed as

```math
    \\|x - y\\|_p = \\left( \\sum_{i=1}^n | x_i - y_i |^p \\right)^{1/p}
```

### Input

- `x` -- vector
- `y` -- vector
- `p` -- (optional, default: `2.0`) the `p`-norm used; `p = 2.0` corresponds to
         the usual Euclidean norm

### Output

A scalar representing ``‖ x - y ‖_p``.
"""
function distance(x::AbstractVector, y::AbstractVector; p::Real=2.0)
    return norm(x - y, p)
end

"""
    allequal(x)

Check whether all elements in a sequence are equal

### Input

- `x` -- sequence

### Output

`true` iff all elements in `x` are equal.

### Notes

The code is taken from [here](https://stackoverflow.com/a/47578613).
"""
function allequal(x)
    length(x) < 2 && return true
    e1 = @inbounds x[1]
    i = 2
    @inbounds for i=2:length(x)
        x[i] == e1 || return false
    end
    return true
end

# if `vector` has exactly one non-zero entry, return its index
# otherwise return 0
function find_unique_nonzero_entry(vector::AbstractVector{N}) where {N}
    res = 0
    for (i, v) in enumerate(vector)
        if v != zero(N)
            if res != 0
                # at least two non-zero entries
                return 0
            else
                # first non-zero entry so far
                res = i
            end
        end
    end
    return res
end

function append_zeros(v::AbstractVector{N}, n::Int) where {N}
    return vcat(v, zeros(N, n))
end

function append_zeros(v::SparseVector{N}, n::Int) where {N}
    return sparsevec(v.nzind, v.nzval, v.n + n)
end

function prepend_zeros(v::AbstractVector{N}, n::Int) where {N}
    return vcat(zeros(N, n), v)
end

function prepend_zeros(v::SparseVector{N}, n::Int) where {N}
    return sparsevec(v.nzind .+ n, v.nzval, v.n + n)
end

"""
    ispermutation(u::AbstractVector{T}, v::AbstractVector) where {T}

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

### Notes

Containment check is performed using `LazySets._in(e, v)`, so in the case of
floating point numbers, the precision to which the check is made is determined
by the type of elements in `v`. See `_in` and `_isapprox` for more information.

Note that approximate equality is not an equivalence relation.
Hence the result may depend on the order of the elements.
"""
function ispermutation(u::AbstractVector{T}, v::AbstractVector) where {T}
    if length(u) != length(v)
        return false
    end
    occurrence_map = Dict{T, Int}()
    has_duplicates = false
    for e in u
        if !_in(e, v)
            return false
        end
        found = false
        for k in keys(occurrence_map)
            if _isapprox(k, e)
                occurrence_map[k] += 1
                has_duplicates = true
                found = true
                break
            end
        end
        if !found
            occurrence_map[e] = 1
        end
    end
    if has_duplicates
        for e in v
            found = false
            for k in keys(occurrence_map)
                if _isapprox(k, e)
                    found = true
                    occurrence_map[k] -= 1
                    if occurrence_map[k] < 0
                        return false
                    end
                    break
                end
            end
            if !found
                return false
            end
        end
    end
    return true
end

"""
    substitute(substitution::Dict{Int, T}, x::AbstractVector{T}) where {T}

Apply a substitution to a given vector.

### Input

- `substitution` -- substitution (a mapping from an index to a new value)
- `x`            -- vector

### Output

A fresh vector corresponding to `x` after `substitution` was applied.
"""
function substitute(substitution::Dict{Int, T}, x::AbstractVector{T}) where {T}
    return substitute!(substitution, copy(x))
end

"""
    substitute!(substitution::Dict{Int, T}, x::AbstractVector{T}) where {T}

Apply a substitution to a given vector.

### Input

- `substitution` -- substitution (a mapping from an index to a new value)
- `x`            -- vector (modified in this function)

### Output

The same (but see the Notes below) vector `x` but after `substitution` was
applied.

### Notes

The vector `x` is modified in-place if it has type `Vector` or `SparseVector`.
Otherwise, we first create a new `Vector` from it.
"""
function substitute!(substitution::Dict{Int, T}, x::AbstractVector{T}) where {T}
    return substitute!(Vector(x), substitution)
end

function substitute!(substitution::Dict{Int, T},
                     x::Union{Vector{T}, SparseVector{T}}) where {T}
    for (index, value) in substitution
        x[index] = value
    end
    return x
end

function _isupwards(vec)
    return vec[2] > 0 || (vec[2] == 0 && vec[1] > 0)
end
