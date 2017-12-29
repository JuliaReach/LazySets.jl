import Base: ∈, ⊆

export AbstractSingleton,
       element,
       an_element

"""
    AbstractSingleton{N<:Real} <: AbstractHyperrectangle{N}

Abstract type for sets with a single value.

### Notes

Every concrete `AbstractSingleton` must define the following functions:
- `element(::AbstractSingleton{N})::Vector{N}` -- return the single element
- `element(::AbstractSingleton{N}, i::Int)::N` -- return the single element's
    entry in the `i`-th dimension

```jldoctest
julia> subtypes(AbstractSingleton)
2-element Array{Union{DataType, UnionAll},1}:
 LazySets.Singleton
 LazySets.ZeroSet
```
"""
abstract type AbstractSingleton{N<:Real} <: AbstractHyperrectangle{N} end


# --- AbstractHyperrectangle interface functions ---


"""
    radius_hyperrectangle(S::AbstractSingleton{N}, i::Int)::N where {N<:Real}

Return the box radius of a set with a single value in a given dimension.

### Input

- `S` -- set with a single value

### Output

Zero.
"""
function radius_hyperrectangle(S::AbstractSingleton{N}, i::Int
                              )::N where {N<:Real}
    return zero(N)
end


"""
    radius_hyperrectangle(S::AbstractSingleton{N})::Vector{N} where {N<:Real}

Return the box radius of a set with a single value in every dimension.

### Input

- `S` -- set with a single value

### Output

The zero vector.
"""
function radius_hyperrectangle(S::AbstractSingleton{N}
                              )::Vector{N} where {N<:Real}
    return zeros(N, dim(S))
end


# --- AbstractPointSymmetric interface functions ---


"""
    center(S::AbstractSingleton{N})::Vector{N} where {N<:Real}

Return the center of a set with a single value.

### Input

- `S` -- set with a single value

### Output

The only element of the set.
"""
function center(S::AbstractSingleton{N})::Vector{N} where {N<:Real}
    return element(S)
end


# --- AbstractPolytope interface functions ---


"""
    vertices_list(S::AbstractSingleton{N})::Vector{Vector{N}} where {N<:Real}

Return the list of vertices of a set with a single value.

### Input

- `S` -- set with a single value

### Output

A list containing only a single vertex.
"""
function vertices_list(S::AbstractSingleton{N}
                      )::Vector{Vector{N}} where {N<:Real}
    return [element(S)]
end


# --- LazySet interface functions ---


"""
    σ(d::AbstractVector{N}, S::AbstractSingleton{N})::Vector{N} where {N<:Real}

Return the support vector of a set with a single value.

### Input

- `d` -- direction
- `S` -- set with a single value

### Output

The support vector, which is the set's vector itself, irrespective of the given
direction.
"""
function σ(d::AbstractVector{N},
           S::AbstractSingleton{N})::Vector{N} where {N<:Real}
    return element(S)
end


"""
    an_element(S::AbstractSingleton{N})::Vector{N} where {N<:Real}

Return some element of a set with a single value.

### Input

- `S` -- set with a single value

### Output

The only element in the set.
"""
function an_element(S::AbstractSingleton{N})::Vector{N} where {N<:Real}
    return element(S)
end


"""
    ∈(x::AbstractVector{N}, S::AbstractSingleton{N})::Bool where {N<:Real}

Check whether a given point is contained in a set with a single value.

### Input

- `x` -- point/vector
- `S` -- set with a single value

### Output

`true` iff ``x ∈ S``.

### Notes

This implementation performs an exact comparison, which may be insufficient with
floating point computations.
"""
function ∈(x::AbstractVector{N}, S::AbstractSingleton{N})::Bool where {N<:Real}
    return x == element(S)
end


"""
    ⊆(S::AbstractSingleton{N}, set::LazySet, witness::Bool=false
     )::Union{Bool,Tuple{Bool,Vector{N}}} where {N<:Real}

Check whether a given set with a single value is contained in a convex set, and
if not, optionally compute a witness.

### Input

- `S`   -- inner set with a single value
- `set` -- outer convex set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``S ⊆ \\text{set}``
* If `witness` option is activated:
  * `(true, [])` iff ``S ⊆ \\text{set}``
  * `(false, v)` iff ``S \\not\\subseteq \\text{set}`` and
    ``v ∈ S \\setminus \\text{set}``
"""
function ⊆(S::AbstractSingleton{N}, set::LazySet, witness::Bool=false
          )::Union{Bool,Tuple{Bool,Vector{N}}} where {N<:Real}
    result = ∈(element(S), set)
    if witness
        return (result, result ? N[] : element(S))
    else
        return result
    end
end


"""
    ⊆(S1::AbstractSingleton{N}, S2::AbstractSingleton{N}, witness::Bool=false
     )::Union{Bool,Tuple{Bool,Vector{N}}} where {N<:Real}

Check whether a given set with a single value is contained in another set with a
single value, and if not, optionally compute a witness.

### Input

- `S1` -- inner set with a single value
- `S2` -- outer set with a single value
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

`true` .

* If `witness` option is deactivated: `true` iff ``S1 ⊆ S2`` iff ``S1 == S2``
* If `witness` option is activated:
  * `(true, [])` iff ``S1 ⊆ S2``
  * `(false, v)` iff ``S1 \\not\\subseteq S2`` and ``v ∈ S1 \\setminus S2``
"""
function ⊆(S1::AbstractSingleton{N},
           S2::AbstractSingleton{N},
           witness::Bool=false
          )::Union{Bool,Tuple{Bool,Vector{N}}} where {N<:Real}
    result = element(S1) == element(S2)
    if witness
        return (result, result ? N[] : element(S1))
    else
        return result
    end
end
