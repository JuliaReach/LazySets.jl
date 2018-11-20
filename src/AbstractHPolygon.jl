import Base: rand,
             ∈

export AbstractHPolygon,
       an_element,
       isredundant,
       remove_redundant_constraints!,
       addconstraint!,
       vertices_list,
       constraints_list,
       isbounded

# This constant marks the threshold for the number of constraints of a polygon
# above which we use a binary search to find the relevant constraint in a
# support vector query.
#
# NOTE: The value must be strictly greater than 2.
const BINARY_SEARCH_THRESHOLD = 10

"""
    AbstractHPolygon{N<:Real} <: AbstractPolygon{N}

Abstract type for polygons in H-representation (i.e., constraints).

### Notes

Every concrete `AbstractHPolygon` must have the following fields:
- `constraints::Vector{LinearConstraint{N}}` -- the constraints

New subtypes should be added to the `convert` method in order to be convertible.

```jldoctest
julia> subtypes(AbstractHPolygon)
2-element Array{Any,1}:
 HPolygon
 HPolygonOpt
```
"""
abstract type AbstractHPolygon{N<:Real} <: AbstractPolygon{N} end


# --- AbstractPolygon interface functions ---


"""
    tovrep(P::AbstractHPolygon{N})::VPolygon{N} where {N<:Real}

Build a vertex representation of the given polygon.

### Input

- `P` -- polygon in constraint representation

### Output

The same polygon but in vertex representation, a `VPolygon`.
"""
function tovrep(P::AbstractHPolygon{N})::VPolygon{N} where {N<:Real}
    return VPolygon(vertices_list(P))
end


"""
    tohrep(P::HPOLYGON)::HPOLYGON where {HPOLYGON<:AbstractHPolygon}

Build a contraint representation of the given polygon.

### Input

- `P` -- polygon in constraint representation

### Output

The identity, i.e., the same polygon instance.
"""
function tohrep(P::HPOLYGON)::HPOLYGON where {HPOLYGON<:AbstractHPolygon}
    return P
end


# --- AbstractPolytope interface functions ---


"""
    vertices_list(P::AbstractHPolygon{N},
                  apply_convex_hull::Bool=false,
                  check_feasibility::Bool=true
                 )::Vector{Vector{N}} where {N<:Real}

Return the list of vertices of a polygon in constraint representation.

### Input

- `P`                 -- polygon in constraint representation
- `apply_convex_hull` -- (optional, default: `false`) flag to post-process the
                         intersection of constraints with a convex hull
- `check_feasibility` -- (optional, default: `true`) flag to check whether the
                         polygon was empty (required for correctness in case of
                         empty polygons)

### Output

List of vertices.

### Algorithm

We compute each vertex as the intersection of consecutive lines defined by the
half-spaces.
If `check_feasibility` is active, we then check if the constraints of the
polygon were actually feasible (i.e., they pointed in the *right* direction).
For this we compute the *average* of all vertices and check membership in each
constraint.
"""
function vertices_list(P::AbstractHPolygon{N},
                       apply_convex_hull::Bool=false,
                       check_feasibility::Bool=true
                      )::Vector{Vector{N}} where {N<:Real}
    n = length(P.constraints)
    points = Vector{Vector{N}}(undef, n)
    if n == 0
        return points
    end
    @inbounds for i in 1:n-1
        points[i] = element(intersection(Line(P.constraints[i]),
                                         Line(P.constraints[i+1])))
    end
    points[n] = element(intersection(Line(P.constraints[n]),
                                     Line(P.constraints[1])))

    # check if polygon was empty
    if check_feasibility
        avg = sum(points) / length(points)
        if avg ∉ P
            return Vector{Vector{N}}(undef, 0)
        end
    end

    return apply_convex_hull ? convex_hull(points) : points
end

"""
    constraints_list(P::AbstractHPolygon{N})::Vector{LinearConstraint{N}} where {N<:Real}

Return the list of constraints defining a polygon in H-representation.

### Input

- `P` -- polygon in H-representation

### Output

The list of constraints of the polygon.
"""
function constraints_list(P::AbstractHPolygon{N})::Vector{LinearConstraint{N}} where {N<:Real}
    return P.constraints
end


# --- LazySet interface functions ---


"""
    an_element(P::AbstractHPolygon{N})::Vector{N} where {N<:Real}

Return some element of a polygon in constraint representation.

### Input

- `P` -- polygon in constraint representation

### Output

A vertex of the polygon in constraint representation (the first one in the order
of the constraints).
"""
function an_element(P::AbstractHPolygon{N})::Vector{N} where {N<:Real}
    @assert length(P.constraints) >= 2 "polygon has less than two constraints"
    return element(intersection(Line(P.constraints[1]),
                                Line(P.constraints[2])))
end

"""
    ∈(x::AbstractVector{N}, P::AbstractHPolygon{N})::Bool where {N<:Real}

Check whether a given 2D point is contained in a polygon in constraint
representation.

### Input

- `x` -- two-dimensional point/vector
- `P` -- polygon in constraint representation

### Output

`true` iff ``x ∈ P``.

### Algorithm

This implementation checks if the point lies on the outside of each edge.
"""
function ∈(x::AbstractVector{N}, P::AbstractHPolygon{N})::Bool where {N<:Real}
    @assert length(x) == 2

    for c in P.constraints
        if dot(c.a, x) > c.b
            return false
        end
    end
    return true
end

"""
    rand(::Type{HPOLYGON}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing,
         [num_constraints]::Int=-1
        )::HPOLYGON{N} where {HPOLYGON<:AbstractHPolygon}

Create a random polygon in constraint representation.

### Input

- `HPOLYGON`        -- type for dispatch
- `N`               -- (optional, default: `Float64`) numeric type
- `dim`             -- (optional, default: 2) dimension
- `rng`             -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`            -- (optional, default: `nothing`) seed for reseeding
- `num_constraints` -- (optional, default: `-1`) number of constraints of the
                       polygon (must be 3 or bigger; see comment below)

### Output

A random polygon in constraint representation.

### Algorithm

We create a random polygon in vertex representation and convert it to constraint
representation.
See [`rand(::Type{VPolygon})`](@ref).
For non-flat polygons the number of vertices and the number of constraints are
identical.
"""
function rand(::Type{HPOLYGON};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int, Nothing}=nothing,
              num_constraints::Int=-1
             )::HPOLYGON{N} where {HPOLYGON<:AbstractHPolygon}
    @assert dim == 2 "cannot create a random $HPOLYGON of dimension $dim"
    @assert num_constraints < 0 || num_constraints >= 3 "cannot construct a " *
        "random $HPOLYGON with only $num_constraints constraints"
    rng = reseed(rng, seed)
    vpolygon = rand(VPolygon; N=N, dim=dim, rng=rng, seed=seed,
                    num_vertices=num_constraints)
    return convert(HPOLYGON, vpolygon)
end


# --- common AbstractHPolygon functions ---

"""
    isredundant(c::LinearConstraint{N}, c1::LinearConstraint{N},
                c2::LinearConstraint{N})::Bool where {N<:Real}

Check whether a linear constraint is redundant wrt. two surrounding constraints.

### Input

- `c`  -- linear constraint of concern
- `c1` -- linear constraint to the right (clockwise turn)
- `c2` -- linear constraint to the left (counter-clockwise turn)

### Output

`true` iff the constraint is redundant.

### Algorithm

We first check whether the angle between the surrounding constraints is < 180°,
which is a necessary condition (unless the direction is identical to one of the
other two constraints).
If so, we next check if the angle is 0°, in which case the constraint `c` is
redundant unless it is strictly tighter than the other two constraints.
If the angle is strictly between 0° and 180°, the constraint `c` is redundant if
and only if the vertex defined by the other two constraints lies inside the set
defined by `c`.
"""
function isredundant(c::LinearConstraint{N},
                     c1::LinearConstraint{N},
                     c2::LinearConstraint{N})::Bool where {N<:Real}
    samedir_check = false
    # determine angle between surrounding constraints
    if !is_right_turn(c1.a, c2.a)
        # angle is > 180°
        samedir_check = true
    elseif is_right_turn(c2.a, c1.a)
        # angle is 0° or 180°
        if samedir(c1.a, c2.a)[1]
            # angle is 0°, i.e., all three constraints have the same direction
            # constraint is redundant unless it is tighter than the other two
            @assert samedir(c1.a, c.a)[1] && samedir(c2.a, c.a)[1]
            return !is_tighter_same_dir_2D(c, c1, strict=true) &&
                   !is_tighter_same_dir_2D(c, c2, strict=true)
        else
            # angle is 180°
            samedir_check = true
        end
    end
    # check if the constraint has the same direction as one of the two
    if samedir(c1.a, c.a)[1]
        return !is_tighter_same_dir_2D(c, c1, strict=true)
    elseif samedir(c2.a, c.a)[1]
        return !is_tighter_same_dir_2D(c, c2, strict=true)
    elseif samedir_check
        # not the same direction => constraint is not redundant
        return false
    end
    cap = intersection(Line(c1), Line(c2))
    @assert cap isa Singleton
    return cap ⊆ c
end

"""
    remove_redundant_constraints!(P::AbstractHPolygon)

Remove all redundant constraints of a polygon in constraint representation.

### Input

- `P` -- polygon in constraint representation

### Output

The same polygon with all redundant constraints removed.

### Notes

Since we only consider bounded polygons and a polygon needs at least three
constraints to be bounded, we stop removing redundant constraints if there are
three or less constraints left.
This means that for non-bounded polygons the result may be unexpected.

### Algorithm

We go through all consecutive triples of constraints and check if the one in the
middle is redundant.
For this we assume that the constraints are sorted.
"""
function remove_redundant_constraints!(P::AbstractHPolygon)
    C = P.constraints
    i = 1
    go_on = true
    c1 = C[length(C)] # define initial c1 here
    while length(C) >= 3 && go_on
        c2 = C[i]
        if i < length(C)
            c3 = C[i+1]
        elseif i == length(C)
            c3 = C[1]
            go_on = false
        end
        if isredundant(c2, c1, c3)
            deleteat!(C, i)
        else
            i += 1
            c1 = C[i-1] # update c1 (updates less often than c2/c3)
        end
    end
    return P
end

"""
    addconstraint!(P::AbstractHPolygon{N},
                   constraint::LinearConstraint{N};
                   [linear_search]::Bool=(length(P.constraints) <
                                          BINARY_SEARCH_THRESHOLD),
                   [prune]::Bool=true
                  )::Nothing where {N<:Real}

Add a linear constraint to a polygon in constraint representation, keeping the
constraints sorted by their normal directions.

### Input

- `P`          -- polygon in constraint representation
- `constraint` -- linear constraint to add
- `linear_search`  -- (optional, default: `length(constraints) <
                      BINARY_SEARCH_THRESHOLD`) flag to choose between linear
                      and binary search
- `prune`          -- (optional, default: `true`) flag for removing redundant
                      constraints in the end

### Output

Nothing.
"""
function addconstraint!(P::AbstractHPolygon{N},
                        constraint::LinearConstraint{N};
                        linear_search::Bool=(length(P.constraints) <
                                             BINARY_SEARCH_THRESHOLD),
                        prune::Bool=true
                       )::Nothing where {N<:Real}
    return addconstraint!(P.constraints, constraint,
                          linear_search=linear_search, prune=prune)
end

"""
    addconstraint!(constraints::Vector{LinearConstraint{N}},
                   new_constraint::LinearConstraint{N};
                   [linear_search]::Bool=(length(P.constraints) <
                                          BINARY_SEARCH_THRESHOLD),
                   [prune]::Bool=true
                  )::Nothing where {N<:Real}

Add a linear constraint to a sorted vector of constrains, keeping the
constraints sorted by their normal directions.

### Input

- `constraints`    -- vector of linear constraintspolygon in constraint
                      representation
- `new_constraint` -- linear constraint to add
- `linear_search`  -- (optional, default: `length(constraints) <
                      BINARY_SEARCH_THRESHOLD`) flag to choose between linear
                      and binary search
- `prune`          -- (optional, default: `true`) flag for removing redundant
                      constraints in the end

### Output

Nothing.

### Algorithm

If `prune` is active, we check if the new constraint is redundant.
If the constraint is not redundant, we perform the same check to the left and to
the right until we find the first constraint that is not redundant.
"""
function addconstraint!(constraints::Vector{LinearConstraint{N}},
                        new_constraint::LinearConstraint{N};
                        linear_search::Bool=(length(constraints) <
                                             BINARY_SEARCH_THRESHOLD),
                        prune::Bool=true
                       )::Nothing where {N<:Real}
    m = length(constraints)
    k = m
    if k > 0
        d = new_constraint.a
        if d <= constraints[1].a
            k = 0
        elseif linear_search
            # linear search
            while d <= constraints[k].a
                k -= 1
            end
        else
            # binary search
            k = binary_search_constraints(
                d, constraints, k, 1 + div(k, 2), choose_lower=true)
        end
    end

    # here constraints[k] <= new_constraint <= constraints[(k%m)+1]
    if prune && m >= 2
        # check if new constraint is redundant
        k += 1
        below = k == 1 ? m : k - 1
        above = k == m + 1 ? 1 : k
        if isredundant(new_constraint, constraints[below], constraints[above])
            return nothing
        end
        # insert new constraint
        insert!(constraints, k, new_constraint)
        m += 1
        # check if old constraints below became redundant
        while m > 2
            center = k == 1 ? m : k - 1
            below = center == 1 ? m : center - 1
            if isredundant(constraints[center], constraints[below],
                           new_constraint)
                deleteat!(constraints, center)
                if center < k
                    k -= 1
                end
                m -= 1
            else
                break
            end
        end
        # check if old constraints above became redundant
        while m > 2
            center = k == m ? 1 : k + 1
            above = center == m ? 1 : center + 1
            if isredundant(constraints[center], new_constraint,
                            constraints[above])
                deleteat!(constraints, center)
                if center < k
                    k -= 1
                end
                m -= 1
            else
                break
            end
        end
    else
        insert!(constraints, k+1, new_constraint)
    end
    return nothing
end

"""
    binary_search_constraints(d::AbstractVector{N},
                              constraints::Vector{LinearConstraint{N}},
                              n::Int,
                              k::Int;
                              [choose_lower]::Bool=false
                             )::Int where {N<:Real}

Performs a binary search in the constraints.

### Input

- `d`            -- direction
- `constraints`  -- constraints
- `n`            -- number of constraints
- `k`            -- start index
- `choose_lower` -- (optional, default: `false`) flag for choosing the lower
                    index (see the 'Output' section)

### Output

In the default setting, the result is the smallest index `k` such that
`d <= constraints[k]`, or `n+1` if no such `k` exists.
If the `choose_lower` flag is set, the result is the largest index `k` such
that `constraints[k] < d`, which is equivalent to being `k-1` in the normal
setting.
"""
function binary_search_constraints(d::AbstractVector{N},
                                   constraints::Vector{LinearConstraint{N}},
                                   n::Int,
                                   k::Int;
                                   choose_lower::Bool=false
                                  )::Int where {N}
    lower = 1
    upper = n+1
    while lower + 1 < upper
        if constraints[k].a <= d
            lower = k
        else
            upper = k
        end
        k = lower + div(upper - lower, 2)
    end
    if choose_lower
        return lower
    else
        if lower == 1 && !(constraints[1].a <= d)
            # special case for index 1
            return 1
        end
        return upper
    end
end

"""
    isbounded(P::AbstractHPolygon, [use_type_assumption]::Bool=true)::Bool

Determine whether a polygon in constraint representation is bounded.

### Input

- `P`                   -- polygon in constraint representation
- `use_type_assumption` -- (optional, default: `true`) flag for ignoring the
                           type assumption that polygons are bounded

### Output

`true` if `use_type_assumption` is activated.
Otherwise, `true` iff `P` is bounded.

### Algorithm

If `!use_type_assumption`, we convert `P` to an `HPolyhedron` `P2` and then use
`isbounded(P2)`.
"""
function isbounded(P::AbstractHPolygon, use_type_assumption::Bool=true)::Bool
    if use_type_assumption
        return true
    end
    return isbounded(HPolyhedron(P.constraints))
end
