import Base.∈

export HalfSpace,
       an_element

"""
    HalfSpace{N<:Real} <: LazySet{N}

Type that represents a (closed) half-space of the form ``a⋅x ≤ b``.

### Fields

- `a` -- normal direction
- `b` -- constraint

### Examples

The set ``y ≥ 0`` in the plane:

```jldoctest
julia> HalfSpace([0, -1.], 0.)
LazySets.HalfSpace{Float64}([0.0, -1.0], 0.0)
```
"""
struct HalfSpace{N<:Real} <: LazySet{N}
    a::Vector{N}
    b::N
end


# --- LazySet interface functions ---


"""
    dim(hs::HalfSpace)::Int

Return the dimension of a half-space.

### Input

- `hs` -- half-space

### Output

The ambient dimension of the half-space.
"""
function dim(hs::HalfSpace)::Int
    return length(hs.a)
end

"""
    σ(d::AbstractVector{<:Real}, hs::HalfSpace)::AbstractVector{<:Real}

Return the support vector of a half-space.

### Input

- `d`  -- direction
- `hs` -- half-space

### Output

The support vector in the given direction, which is only defined in the
following two cases:
1. The direction has norm zero.
2. The direction is the half-space's normal direction.
In both cases the result is any point on the boundary (the defining hyperplane).
Otherwise this function throws an error.
"""
function σ(d::AbstractVector{N},
           hs::HalfSpace)::AbstractVector{<:Real} where {N<:Real}
    @assert (length(d) == dim(hs)) "cannot compute the support vector of a " *
        "$(dim(hs))-dimensional half-space along a vector of length $(length(d))"

    @inline function not_solvable(d, hs)
        error("the support vector for the half-space with normal direction " *
            "$(hs.a) is not defined along a direction $d")
    end

    first_nonzero_entry_a = -1
    if all(d .== 0)
        # zero vector
        return an_element(hs)
    else
        # not the zero vector, check if it is a normal vector
        factor = zero(N)
        for i in 1:length(hs.a)
            if hs.a[i] == 0
                if d[i] != 0
                    not_solvable(d, hs)
                end
            else
                if d[i] == 0
                    not_solvable(d, hs)
                elseif first_nonzero_entry_a == -1
                    factor = hs.a[i] / d[i]
                    first_nonzero_entry_a = i
                elseif d[i] * factor != hs.a[i]
                    not_solvable(d, hs)
                end
            end
        end
        return an_element_helper(hs, first_nonzero_entry_a)
    end
end

"""
    an_element(hs::HalfSpace{N})::Vector{N} where {N<:Real}

Return some element of a half-space.

### Input

- `hs` -- half-space

### Output

An element in the half-space.
"""
function an_element(hs::HalfSpace{N})::Vector{N} where {N<:Real}
    return an_element_helper(hs, findfirst(hs.a))
end

"""
    ∈(x::AbstractVector{N}, hs::HalfSpace{N})::Bool where {N<:Real}

Check whether a given point is contained in a half-space.

### Input

- `x` -- point/vector
- `hs` -- half-space

### Output

`true` iff ``x ∈ hs``.

### Algorithm

We just check if ``x`` satisfies ``a⋅x ≤ b``.
"""
function ∈(x::AbstractVector{N}, hs::HalfSpace{N})::Bool where {N<:Real}
    return dot(x, hs.a) <= hs.b
end


# --- HalfSpace functions ---


"""
    an_element_helper(hs::HalfSpace{N},
                      first_nonzero_entry_a::Int)::Vector{N} where {N<:Real}

Helper function that computes an element on a half-space's hyperplane.

### Input

- `hs` -- half-space
- `first_nonzero_entry_a` -- index such that `hs.a` is different from 0

### Output

An element on a half-space's hyperplane.

### Algorithm

We compute the point on the hyperplane as follows:
- We already found a nonzero entry of ``a`` in dimension, say, ``i``.
- We set ``x[j] = 0`` for ``j ≠ i``.
- We set ``x[i] = b / a[i]``.
"""
function an_element_helper(hs::HalfSpace{N},
                           first_nonzero_entry_a::Int)::Vector{N} where {N<:Real}
    @assert first_nonzero_entry_a in 1:length(hs.a) "invalid index " *
        "$first_nonzero_entry_a for half-space"
    x = zeros(N, dim(hs))
    x[first_nonzero_entry_a] = hs.b / hs.a[first_nonzero_entry_a]
    return x
end
