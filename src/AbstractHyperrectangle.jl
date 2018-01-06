import Base.LinAlg.norm,
       Base.∈,
       Base.⊆,
       Base.∩

export AbstractHyperrectangle,
       radius_hyperrectangle

"""
    AbstractHyperrectangle{N<:Real} <: AbstractPointSymmetricPolytope{N}

Abstract type for box-shaped sets.

### Notes

Every concrete `AbstractHyperrectangle` must define the following functions:
- `radius_hyperrectangle(::AbstractHyperrectangle{N})::Vector{N}` -- return the
    hyperrectangle's radius, which is a full-dimensional vector
- `radius_hyperrectangle(::AbstractHyperrectangle{N}, i::Int)::N` -- return the
    hyperrectangle's radius in the `i`-th dimension

```jldoctest
julia> subtypes(AbstractHyperrectangle)
3-element Array{Union{DataType, UnionAll},1}:
 LazySets.AbstractSingleton
 LazySets.BallInf       
 LazySets.Hyperrectangle
```
"""
abstract type AbstractHyperrectangle{N<:Real} <: AbstractPointSymmetricPolytope{N}
end


# --- AbstractPolytope interface functions ---


"""
    vertices_list(B::AbstractHyperrectangle{N})::Vector{Vector{N}} where {N<:Real}

Return the list of vertices of a box-shaped set.

### Input

- `B` -- box-shaped set

### Output

A list of vertices.

### Notes

For high dimensions, it is preferable to develop a `vertex_iterator` approach.
"""
function vertices_list(B::AbstractHyperrectangle{N}
                      )::Vector{Vector{N}} where {N<:Real}
    return [center(B) .+ si .* radius_hyperrectangle(B)
        for si in IterTools.product([[1, -1] for i = 1:dim(B)]...)]
end


# --- LazySet interface functions ---


"""
    σ(d::AbstractVector{N}, B::AbstractHyperrectangle{N}
     )::AbstractVector{N} where {N<:Real}

Return the support vector of a box-shaped set in a given direction.

### Input

- `d` -- direction
- `B` -- box-shaped set

### Output

The support vector in the given direction.
If the direction has norm zero, the vertex with biggest values is returned.
"""
function σ(d::AbstractVector{N},
           B::AbstractHyperrectangle{N})::AbstractVector{N} where {N<:Real}
    return center(B) .+ sign_cadlag.(d) .* radius_hyperrectangle(B)
end

"""
    norm(B::AbstractHyperrectangle, [p]::Real=Inf)::Real

Return the norm of a box-shaped set.

### Input

- `B` -- box-shaped set
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the norm.

### Notes

The norm of a box-shaped set is defined as the norm of the enclosing ball, of
the given ``p``-norm, of minimal volume.
"""
function norm(B::AbstractHyperrectangle, p::Real=Inf)::Real
    return maximum(map(x -> norm(x, p), vertices_list(B)))
end

"""
    diameter(B::AbstractHyperrectangle, [p]::Real=Inf)::Real

Return the diameter of a box-shaped set.

### Input

- `H` -- box-shaped set
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the diameter.

### Notes

The diameter is defined as the maximum distance in the given ``p``-norm between
any two elements of the set.
Equivalently, it is the diameter of the enclosing ball of the given ``p``-norm
of minimal volume with the same center.
"""
function diameter(B::AbstractHyperrectangle, p::Real=Inf)::Real
    return radius(B, p) * 2
end

"""
    ∈(x::AbstractVector{N}, B::AbstractHyperrectangle{N})::Bool where {N<:Real}

Check whether a given point is contained in a box-shaped set.

### Input

- `x` -- point/vector
- `B` -- box-shaped set

### Output

`true` iff ``x ∈ B``.

### Algorithm

Let ``B`` be an ``n``-dimensional box-shaped set, ``c_i`` and ``r_i`` be
the box's center and radius and ``x_i`` be the vector ``x`` in dimension ``i``,
respectively.
Then ``x ∈ B`` iff ``|c_i - x_i| ≤ r_i`` for all ``i=1,…,n``.
"""
function ∈(x::AbstractVector{N},
           B::AbstractHyperrectangle{N})::Bool where {N<:Real}
    @assert length(x) == dim(B)
    for i in eachindex(x)
        if abs(center(B)[i] - x[i]) > radius_hyperrectangle(B, i)
            return false
        end
    end
    return true
end

"""
    ⊆(S::LazySet, H::AbstractHyperrectangle{N}, witness::Bool=false
     )::Union{Bool,Tuple{Bool,Vector{N}}} where {N<:Real}

Check whether a convex set is contained in a hyperrectangle, and if not,
optionally compute a witness.

### Input

- `S` -- inner convex set
- `H` -- outer hyperrectangle
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``S ⊆ H``
* If `witness` option is activated:
  * `(true, [])` iff ``S ⊆ H``
  * `(false, v)` iff ``S \\not\\subseteq H`` and ``v ∈ S \\setminus H``

### Algorithm

``S ⊆ H`` iff ``\\operatorname{ihull}(S) ⊆ H``, where  ``\\operatorname{ihull}``
is the interval hull operator.
"""
function ⊆(S::LazySet, H::AbstractHyperrectangle{N}, witness::Bool=false
          )::Union{Bool,Tuple{Bool,Vector{N}}} where {N<:Real}
    return ⊆(Approximations.interval_hull(S), H, witness)
end

"""
    ⊆(P::AbstractPolytope{N}, H::AbstractHyperrectangle, witness::Bool=false
     )::Union{Bool,Tuple{Bool,Vector{N}}} where {N<:Real}

Check whether a polytope is contained in a hyperrectangle, and if not,
optionally compute a witness.

### Input

- `P` -- inner polytope
- `H` -- outer hyperrectangle
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``P ⊆ H``
* If `witness` option is activated:
  * `(true, [])` iff ``P ⊆ H``
  * `(false, v)` iff ``P \\not\\subseteq H`` and ``v ∈ P \\setminus H``

### Notes

This copy-pasted method just exists to avoid method ambiguities.

### Algorithm

Since ``H`` is convex, ``P ⊆ H`` iff ``v_i ∈ H`` for all vertices ``v_i`` of
``P``.
"""
function ⊆(P::AbstractPolytope{N},
           H::AbstractHyperrectangle,
           witness::Bool=false
          )::Union{Bool,Tuple{Bool,Vector{N}}} where {N<:Real}
    @assert dim(P) == dim(H)

    for v in vertices_list(P)
        if !∈(v, H)
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

"""
    ⊆(H1::AbstractHyperrectangle{N}, H2::AbstractHyperrectangle{N},
      witness::Bool=false)::Union{Bool,Tuple{Bool,Vector{N}}} where {N<:Real}

Check whether a given hyperrectangle is contained in another hyperrectangle, and
if not, optionally compute a witness.

### Input

- `H1` -- inner hyperrectangle
- `H2` -- outer hyperrectangle
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``H1 ⊆ H2``
* If `witness` option is activated:
  * `(true, [])` iff ``H1 ⊆ H2``
  * `(false, v)` iff ``H1 \\not\\subseteq H2`` and ``v ∈ H1 \\setminus H2``

### Algorithm

``H1 ⊆ H2`` iff ``c_1 + r_1 ≤ c_2 + r_2 ∧ c_1 - r_1 ≥ c_2 - r_2`` iff
``r_1 - r_2 ≤ c_1 - c_2 ≤ -(r_1 - r_2)``, where ``≤`` is taken component-wise.
"""
function ⊆(H1::AbstractHyperrectangle{N},
           H2::AbstractHyperrectangle{N},
           witness::Bool=false
          )::Union{Bool,Tuple{Bool,Vector{N}}} where {N<:Real}
    @assert dim(H1) == dim(H2)

    for i in 1:dim(H1)
        c_dist = center(H1)[i] - center(H2)[i]
        r_dist = radius_hyperrectangle(H1, i) - radius_hyperrectangle(H2, i)
        if -r_dist < c_dist || c_dist < r_dist
            if witness
                # compute a witness 'p' in the difference
                p = copy(center(H1))
                if c_dist >= 0
                    p[i] += radius_hyperrectangle(H1, i)
                else
                    p[i] -= radius_hyperrectangle(H1, i)
                end
                return (false, p)
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

"""
    ∩(H1::AbstractHyperrectangle{N}, H2::AbstractHyperrectangle{N},
      witness::Bool=false)::Union{Bool,Tuple{Bool,Vector{N}}} where {N<:Real}

Check whether two hyperrectangles intersect, and if so, optionally compute a
witness.

### Input

- `H1` -- first hyperrectangle
- `H2` -- second hyperrectangle
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``H1 ∩ H2 ≠ ∅``
* If `witness` option is activated:
  * `(true, v)` iff ``H1 ∩ H2 ≠ ∅`` and ``v ∈ H1 ∩ H2``
  * `(false, [])` iff ``H1 ∩ H2 = ∅``

### Algorithm

``H1 ∩ H2 ≠ ∅`` iff ``|c_2 - c_1| ≤ r_1 + r_2``, where ``≤`` is taken
component-wise.

A witness is computed by starting in one center and moving toward the other
center for as long as the minimum of the radius and the center distance.
In other words, the witness is the point in `H1` that is closest to the center
of `H2`.
"""
function ∩(H1::AbstractHyperrectangle{N},
           H2::AbstractHyperrectangle{N},
           witness::Bool=false
          )::Union{Bool,Tuple{Bool,Vector{N}}} where {N<:Real}
    intersection = true
    center_diff = center(H2) - center(H1)
    for i in eachindex(center_diff)
        if abs(center_diff[i]) >
                radius_hyperrectangle(H1, i) + radius_hyperrectangle(H2, i)
            intersection = false
            break
        end
    end

    if !witness
        return intersection
    elseif !intersection
        return (false, N[])
    end

    # compute a witness 'p' in the intersection
    p = copy(center(H1))
    c2 = center(H2)
    for i in eachindex(center_diff)
        if p[i] <= c2[i]
            # second center is right of first center
            p[i] += min(radius_hyperrectangle(H1, i), center_diff[i])
        else
            # second center is left of first center
            p[i] -= max(radius_hyperrectangle(H1, i), -center_diff[i])
        end
    end
    return (true, p)
end
