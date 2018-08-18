import Base.∈

export Hyperplane,
       an_element

"""
    Hyperplane{N<:Real} <: LazySet{N}

Type that represents a hyperplane of the form ``a⋅x = b``.

### Fields

- `a` -- normal direction
- `b` -- constraint

### Examples

The plane ``y = 0``:

```jldoctest
julia> Hyperplane([0, 1.], 0.)
Hyperplane{Float64}([0.0, 1.0], 0.0)
```
"""
struct Hyperplane{N<:Real} <: LazySet{N}
    a::AbstractVector{N}
    b::N
end


# --- LazySet interface functions ---


"""
    dim(hp::Hyperplane)::Int

Return the dimension of a hyperplane.

### Input

- `hp` -- hyperplane

### Output

The ambient dimension of the hyperplane.
"""
function dim(hp::Hyperplane)::Int
    return length(hp.a)
end

"""
    σ(d::AbstractVector{N}, hp::Hyperplane{N}) where {N<:Real}

Return the support vector of a hyperplane.

### Input

- `d`  -- direction
- `hp` -- hyperplane

### Output

The support vector in the given direction, which is only defined in the
following two cases:
1. The direction has norm zero.
2. The direction is the hyperplane's normal direction.
In both cases the result is any point on the hyperplane.
Otherwise this function throws an error.
"""
function σ(d::AbstractVector{N}, hp::Hyperplane{N}) where {N<:Real}
    return σ_helper(d, hp)
end

"""
    an_element(hp::Hyperplane{N})::Vector{N} where {N<:Real}

Return some element of a hyperplane.

### Input

- `hp` -- hyperplane

### Output

An element on the hyperplane.
"""
function an_element(hp::Hyperplane{N})::Vector{N} where {N<:Real}
    return an_element_helper(hp)
end

"""
    ∈(x::AbstractVector{N}, hp::Hyperplane{N})::Bool where {N<:Real}

Check whether a given point is contained in a hyperplane.

### Input

- `x` -- point/vector
- `hp` -- hyperplane

### Output

`true` iff ``x ∈ hp``.

### Algorithm

We just check if ``x`` satisfies ``a⋅x = b``.
"""
function ∈(x::AbstractVector{N}, hp::Hyperplane{N})::Bool where {N<:Real}
    return dot(x, hp.a) == hp.b
end


# --- Hyperplane functions ---


"""
```
    σ_helper(d::AbstractVector{N},
             hp::Hyperplane{N},
             [name]::String="hyperplane") where {N<:Real}
```

Return the support vector of a hyperplane.

### Input

- `d`  -- direction
- `hp` -- hyperplane
- `name` -- (optional, default: "hyperplane") name for error messages

### Output

The support vector in the given direction, which is only defined in the
following two cases:
1. The direction has norm zero.
2. The direction is the hyperplane's normal direction.
In both cases the result is any point on the hyperplane.
Otherwise this function throws an error.
"""
@inline function σ_helper(d::AbstractVector{N},
                          hp::Hyperplane{N},
                          name::String="hyperplane") where {N<:Real}
    @assert (length(d) == dim(hp)) "cannot compute the support vector of a " *
        "$(dim(hp))-dimensional $name along a vector of length $(length(d))"

    @inline function not_solvable(d, hp)
        error("the support vector for the $name with normal direction " *
            "$(hp.a) is not defined along a direction $d")
    end

    first_nonzero_entry_a = -1
    if all(d .== 0)
        # zero vector
        return an_element(hp)
    else
        # not the zero vector, check if it is a normal vector
        factor = zero(N)
        for i in 1:length(hp.a)
            if hp.a[i] == 0
                if d[i] != 0
                    not_solvable(d, hp)
                end
            else
                if d[i] == 0
                    not_solvable(d, hp)
                elseif first_nonzero_entry_a == -1
                    factor = hp.a[i] / d[i]
                    first_nonzero_entry_a = i
                elseif d[i] * factor != hp.a[i]
                    not_solvable(d, hp)
                end
            end
        end
        return an_element_helper(hp, first_nonzero_entry_a)
    end
end

"""
    an_element_helper(hp::Hyperplane{N},
                      [nonzero_entry_a]::Int)::Vector{N} where {N<:Real}

Helper function that computes an element on a hyperplane's hyperplane.

### Input

- `hp` -- hyperplane
- `nonzero_entry_a` -- (optional, default: computes the first index) index `i`
                       such that `hp.a[i]` is different from 0

### Output

An element on a hyperplane.

### Algorithm

We compute the point on the hyperplane as follows:
- We already found a nonzero entry of ``a`` in dimension, say, ``i``.
- We set ``x[i] = b / a[i]``.
- We set ``x[j] = 0`` for all ``j ≠ i``.
"""
@inline function an_element_helper(hp::Hyperplane{N},
                                   nonzero_entry_a::Int=findnext(x -> x!=zero(N), hp.a, 1)
                                  )::Vector{N} where {N<:Real}
    @assert nonzero_entry_a in 1:length(hp.a) "invalid index " *
        "$nonzero_entry_a for hyperplane"
    x = zeros(N, dim(hp))
    x[nonzero_entry_a] = hp.b / hp.a[nonzero_entry_a]
    return x
end
