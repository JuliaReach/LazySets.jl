import Base.⊆

export AbstractPolytope,
       vertices_list,
       singleton_list

"""
    AbstractPolytope{N<:Real} <: LazySet{N}

Abstract type for polytopic sets, i.e., sets with finitely many flat facets, or
equivalently, sets defined as an intersection of a finite number of halfspaces,
or equivalently, sets with finitely many vertices.

### Notes

Every concrete `AbstractPolytope` must define the following functions:
- `vertices_list(::AbstractPolytope{N})::Vector{Vector{N}}` -- return a list of
    all vertices

```jldoctest
julia> subtypes(AbstractPolytope)
3-element Array{Union{DataType, UnionAll},1}:
 LazySets.AbstractPointSymmetricPolytope
 LazySets.AbstractPolygon
 LazySets.HPolytope
```
"""
abstract type AbstractPolytope{N<:Real} <: LazySet{N} end


# --- LazySet interface functions ---


"""
    ⊆(P::AbstractPolytope{N}, S::LazySet{N}, witness::Bool=false
     )::Union{Bool,Tuple{Bool,Vector{N}}} where {N<:Real}

Check whether a polytope is contained in a convex set, and if not, optionally
compute a witness.

### Input

- `P` -- inner polytope
- `S` -- outer convex set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``P ⊆ S``
* If `witness` option is activated:
  * `(true, [])` iff ``P ⊆ S``
  * `(false, v)` iff ``P \\not\\subseteq S`` and ``v ∈ P \\setminus S``

### Algorithm

Since ``S`` is convex, ``P ⊆ S`` iff ``v_i ∈ S`` for all vertices ``v_i`` of
``P``.
"""
function ⊆(P::AbstractPolytope{N}, S::LazySet{N}, witness::Bool=false
          )::Union{Bool,Tuple{Bool,Vector{N}}} where {N<:Real}
    @assert dim(P) == dim(S)

    for v in vertices_list(P)
        if !∈(v, S)
            if witness
                return (false, v)
            else
                return false
            end
        end
    end
    if witness
        return (true, N[])
    else
        return true
    end
end


# --- common AbstractPolytope functions ---


"""
    singleton_list(P::AbstractPolytope{N})::Vector{Singleton{N}} where {N<:Real}

Return the vertices of a polytopic as a list of singletons.

### Input

- `P` -- a polytopic set

### Output

List containing a singleton for each vertex.
"""
function singleton_list(P::AbstractPolytope{N}
                       )::Vector{Singleton{N}} where {N<:Real}
    return [Singleton(vi) for vi in P.vertices_list]
end
