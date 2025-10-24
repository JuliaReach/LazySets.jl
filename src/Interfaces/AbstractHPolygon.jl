export AbstractHPolygon,
       isredundant

# This constant marks the threshold for the number of constraints of an
# H-polygon above which we use a binary search to find the relevant constraint
# in a support-vector query.
#
# NOTE: The value must be strictly greater than 2.
const BINARY_SEARCH_THRESHOLD = 10

"""
    AbstractHPolygon{N} <: AbstractPolygon{N}

Abstract type for polygons in constraint representation.

### Notes

See [`HPolygon`](@ref) for a standard implementation of this interface.

All subtypes must satisfy the invariant that constraints are sorted
counter-clockwise.

Every concrete `AbstractHPolygon` must have the following fields:

- `constraints::Vector{HalfSpace{N,AbstractVector{N}}}` -- constraints vector

The subtypes of `AbstractHPolygon`:

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractHPolygon)
2-element Vector{Any}:
 HPolygon
 HPolygonOpt
```
"""
abstract type AbstractHPolygon{N} <: AbstractPolygon{N} end

"""
    tovrep(P::AbstractHPolygon)

Build a vertex representation of a polygon in constraint representation.

### Input

- `P` -- polygon in constraint representation

### Output

The same polygon but in vertex representation, a `VPolygon`.
"""
function tovrep(P::AbstractHPolygon)
    return VPolygon(vertices_list(P; apply_convex_hull=false))
end

"""
    tohrep(P::HPOLYGON) where {HPOLYGON<:AbstractHPolygon}

Build a constraint representation of the given polygon.

### Input

- `P` -- polygon in constraint representation

### Output

The identity, i.e., the same polygon instance.
"""
function tohrep(P::HPOLYGON) where {HPOLYGON<:AbstractHPolygon}
    return P
end

"""
    normalize(P::AbstractHPolygon{N}, p::Real=N(2)) where {N}

Normalize a polygon in constraint representation.

### Input

- `P` -- polygon in constraint representation
- `p` -- (optional, default: `2`) norm

### Output

A new polygon in constraint representation whose normal directions ``a_i``
are normalized, i.e., such that ``‖a_i‖_p = 1`` holds.
"""
function normalize(P::AbstractHPolygon{N}, p::Real=N(2)) where {N}
    constraints = [normalize(hs, p) for hs in constraints_list(P)]
    T = basetype(P)
    return T(constraints)
end

"""
    vertices_list(P::AbstractHPolygon;
                  apply_convex_hull::Bool=true,
                  check_feasibility::Bool=true)

Return the list of vertices of a polygon in constraint representation.

### Input

- `P`                 -- polygon in constraint representation
- `apply_convex_hull` -- (optional, default: `true`) flag to post-process the
                         intersection of constraints with a convex hull
- `check_feasibility` -- (optional, default: `true`) flag to check whether the
                         polygon was empty (required for correctness in case of
                         empty polygons)

### Output

List of vertices.

### Notes

By construction an `AbstractHPolygon` should not contain any redundant vertices.
Still the `apply_convex_hull` argument is activated by default to remove
potential duplicate vertices. They can exist due to numeric instability.

```jldoctest
julia> p = HPolygon([HalfSpace([1.0, 0.0], 1.0),
                     HalfSpace([0.0, 1.0], 1.0),
                     HalfSpace([-1.0, 0.0], -1.0),
                     HalfSpace([0.0, -1.0], -1.0)]);

julia> vertices_list(p, apply_convex_hull=false)
4-element Vector{Vector{Float64}}:
 [1.0, 1.0]
 [1.0, 1.0]
 [1.0, 1.0]
 [1.0, 1.0]
```

If it is known that each constraint has a "proper" distance to the next vertex,
this step can be skipped.

### Algorithm

We compute each vertex as the intersection of consecutive lines defined by the
half-spaces.
If `check_feasibility` is active, we then check if the constraints of the
polygon were actually feasible (i.e., they pointed in the *right* direction).
For this we compute the *average* of all vertices and check membership in each
constraint.
"""
function vertices_list(P::AbstractHPolygon;
                       apply_convex_hull::Bool=true,  # see docstring and #1405
                       check_feasibility::Bool=true)
    n = length(P.constraints)
    N = eltype(P)
    points = Vector{Vector{N}}(undef, n)
    if n == 0
        throw(ArgumentError("a polygon with no constraints is invalid"))
    end
    @inbounds for i in 1:(n - 1)
        cap = _intersection_line2d(P.constraints[i], P.constraints[i + 1])
        if cap isa EmptySet
            return Vector{Vector{N}}()
        else
            points[i] = element(cap)
        end
    end
    cap = _intersection_line2d(P.constraints[n], P.constraints[1])
    if cap isa EmptySet
        return Vector{Vector{N}}()
    else
        points[n] = element(cap)
    end

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
    constraints_list(P::AbstractHPolygon)

Return the list of constraints defining a polygon in constraint representation.

### Input

- `P` -- polygon in constraint representation

### Output

The list of constraints of the polygon.
The implementation guarantees that the constraints are sorted counter-clockwise.
"""
@validate function constraints_list(P::AbstractHPolygon)
    return P.constraints
end

"""
    an_element(P::AbstractHPolygon)

Return some element of a polygon in constraint representation.

### Input

- `P` -- polygon in constraint representation

### Output

A vertex of the polygon in constraint representation (the first one in the order
of the constraints).
"""
function an_element(P::AbstractHPolygon)
    @assert length(P.constraints) >= 2 "polygon has less than two constraints"
    return element(_intersection_line2d(P.constraints[1], P.constraints[2]))
end

"""
    ∈(x::AbstractVector, P::AbstractHPolygon)

Check whether a given two-dimensional point is contained in a polygon in
constraint representation.

### Input

- `x` -- two-dimensional point/vector
- `P` -- polygon in constraint representation

### Output

`true` iff ``x ∈ P``.

### Algorithm

This implementation checks if the point lies inside each constraint.
"""
@validate function ∈(x::AbstractVector, P::AbstractHPolygon)
    for c in P.constraints
        if !_leq(dot(c.a, x), c.b)
            return false
        end
    end
    return true
end

"""
    rand(::Type{HPOLYGON}; [N]::Type=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing,
         [num_constraints]::Int=-1) where {HPOLYGON<:AbstractHPolygon}

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
              N::Type=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing,
              num_constraints::Int=-1) where {HPOLYGON<:AbstractHPolygon}
    @assert dim == 2 "cannot create a random $HPOLYGON of dimension $dim"
    @assert num_constraints < 0 || num_constraints >= 3 "cannot construct a " *
                                                        "random $HPOLYGON with only $num_constraints constraints"
    rng = reseed!(rng, seed)
    vpolygon = rand(VPolygon; N=N, dim=dim, rng=rng, seed=seed,
                    num_vertices=num_constraints)
    return convert(HPOLYGON, vpolygon)
end

"""
    isredundant(cmid::HalfSpace, cright::HalfSpace, cleft::HalfSpace)

Check whether a linear constraint is redundant wrt. two surrounding constraints.

### Input

- `cmid`   -- linear constraint of concern
- `cright` -- linear constraint to the right (clockwise turn)
- `cleft`  -- linear constraint to the left (counter-clockwise turn)

### Output

`true` iff the constraint is redundant.

### Algorithm

We first check whether the angle between the surrounding constraints is < 180°,
which is a necessary condition (unless the direction is identical to one of the
other two constraints).
If so, we next check if the angle is 0°, in which case the constraint `cmid` is
redundant unless it is strictly tighter than the other two constraints.
If the angle is strictly between 0° and 180°, the constraint `cmid` is redundant
if and only if the vertex defined by the other two constraints lies inside the
set defined by `cmid`.
"""
function isredundant(cmid::HalfSpace, cright::HalfSpace, cleft::HalfSpace)
    samedir_check = false
    # determine angle between surrounding constraints
    if !is_right_turn(cright.a, cleft.a)
        # angle is > 180°
        samedir_check = true
    elseif is_right_turn(cleft.a, cright.a)
        # angle is 0° or 180°
        if samedir(cright.a, cleft.a)[1]
            # angle is 0°
            if samedir(cright.a, cmid.a)[1] && samedir(cleft.a, cmid.a)[1]
                # all three constraints have the same direction
                # constraint is redundant unless it is tighter than the others
                return !is_tighter_same_dir_2D(cmid, cright; strict=true) &&
                       !is_tighter_same_dir_2D(cmid, cleft; strict=true)
            else
                # surrounding constraints have the same direction but the
                # central constraint does not => corner case with just three
                # non-bounding constraints (usually occurs during incremental
                # addition of constraints)
                return false
            end
        else
            # angle is 180°
            samedir_check = true
        end
    end
    # check if the constraint has the same direction as one of the two
    if samedir(cright.a, cmid.a)[1]
        return !is_tighter_same_dir_2D(cmid, cright; strict=true)
    elseif samedir(cleft.a, cmid.a)[1]
        return !is_tighter_same_dir_2D(cmid, cleft; strict=true)
    elseif samedir_check
        # not the same direction => constraint is not redundant
        return false
    end
    cap = _intersection_line2d(cright, cleft)
    @assert cap isa Singleton
    return cap ⊆ cmid
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
three or fewer constraints left.
Hence for unbounded polygons the result may be unexpected.

### Algorithm

We go through all consecutive triples of constraints and check if the one in the
middle is redundant.
For this we assume that the constraints are sorted.
"""
function remove_redundant_constraints!(P::AbstractHPolygon)
    C = P.constraints
    i = 1
    go_on = true
    cright = C[length(C)] # define initial cright here
    while length(C) >= 3 && go_on
        cmid = C[i]
        if i < length(C)
            cleft = C[i + 1]
        elseif i == length(C)
            cleft = C[1]
            go_on = false
        end
        if isredundant(cmid, cright, cleft)
            deleteat!(C, i)
        else
            i += 1
            cright = C[i - 1] # update cright (updates less often than cmid/cleft)
        end
    end
    return P
end

"""
    addconstraint!(P::AbstractHPolygon, constraint::HalfSpace;
                   [linear_search]::Bool=length(P.constraints) < $BINARY_SEARCH_THRESHOLD,
                   [prune]::Bool=true)

Add a linear constraint to a polygon in constraint representation, keeping the
constraints sorted by their normal directions.

### Input

- `P`             -- polygon in constraint representation
- `constraint`    -- linear constraint to add
- `linear_search` -- (optional, default:
                     `length(constraints) < $BINARY_SEARCH_THRESHOLD`) flag to
                     choose between linear and binary search
- `prune`         -- (optional, default: `true`) flag for removing redundant
                     constraints in the end
"""
function addconstraint!(P::AbstractHPolygon, constraint::HalfSpace;
                        linear_search::Bool=length(P.constraints) < BINARY_SEARCH_THRESHOLD,
                        prune::Bool=true)
    return addconstraint!(P.constraints, constraint;
                          linear_search=linear_search, prune=prune)
end

"""
    addconstraint!(constraints::Vector{<:HalfSpace}, new_constraint::HalfSpace;
                   [linear_search]::Bool=length(P.constraints) < $BINARY_SEARCH_THRESHOLD,
                   [prune]::Bool=true)

Add a linear constraint to a sorted vector of constrains, keeping the
constraints sorted by their normal directions.

### Input

- `constraints`    -- vector of linear constraints
- `new_constraint` -- linear constraint to add
- `linear_search`  -- (optional, default:
                      `length(constraints) < $BINARY_SEARCH_THRESHOLD`) flag to
                      choose between linear and binary search
- `prune`          -- (optional, default: `true`) flag for removing redundant
                      constraints in the end

### Algorithm

If `prune` is active, we check if the new constraint is redundant.
If the constraint is not redundant, we perform the same check to the left and to
the right until we find the first constraint that is not redundant.
"""
function addconstraint!(constraints::Vector{<:HalfSpace}, new_constraint::HalfSpace;
                        linear_search::Bool=length(constraints) < BINARY_SEARCH_THRESHOLD,
                        prune::Bool=true)
    m = length(constraints)
    k = m
    if k > 0
        d = new_constraint.a
        if d ⪯ constraints[1].a
            k = 0
        elseif linear_search
            # linear search
            while d ⪯ constraints[k].a
                k -= 1
            end
        else
            # binary search
            k = binary_search_constraints(d, constraints; choose_lower=true)
        end
    end

    # here constraints[k].a ⪯ new_constraint.a ⪯ constraints[(k%m)+1].a
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
        insert!(constraints, k + 1, new_constraint)
    end
    return nothing
end

"""
    binary_search_constraints(d::AbstractVector{N},
                              constraints::Vector{<:HalfSpace{N}};
                              [start_index]::Int=div(length(constraints)+1, 2),
                              [choose_lower]::Bool=false) where {N}

Perform a binary search in the constraints.

### Input

- `d`            -- direction
- `constraints`  -- constraints
- `start_index`  -- (optional, default: `div(length(constraints)+1, 2)`) start
                    index
- `choose_lower` -- (optional, default: `false`) flag for choosing the lower
                    index (see the 'Output' section)

### Output

In the default setting, the result is the smallest index `k` such that
`d ⪯ constraints[k].a`, or `length(constraints)+1` if no such `k` exists.
If the `choose_lower` flag is set, the result is the largest index `k` such
that `constraints[k].a < d`, which is equivalent to being `k-1` in the normal
setting.
"""
function binary_search_constraints(d::AbstractVector{N},
                                   constraints::Vector{<:HalfSpace{N}};
                                   start_index::Int=div(length(constraints) + 1, 2),
                                   choose_lower::Bool=false) where {N}
    lower = 1
    n = length(constraints)
    upper = n + 1
    @assert 1 <= start_index <= n "invalid start index $start_index"
    m = start_index
    while lower + 1 < upper
        if constraints[m].a ⪯ d
            lower = m
        else
            upper = m
        end
        m = div(lower + upper, 2)
    end

    # since `⪯` is approximate, it can happen that x ⪯ y and y ⪯ x
    # we want to return the smallest index, so (linearly) search to the left
    while lower > 1 && d ⪯ constraints[lower].a && constraints[lower].a ⪯ d
        lower -= 1
        upper -= 1
    end

    if choose_lower
        return lower
    else
        if lower == 1 && !(constraints[1].a ⪯ d)
            # special case: smaller than all elements in the vector
            return 1
        end
        return upper
    end
end

"""
    isbounded(P::AbstractHPolygon, [use_type_assumption]::Bool=true)

Determine whether a polygon in constraint representation is bounded.

### Input

- `P`                   -- polygon in constraint representation
- `use_type_assumption` -- (optional, default: `true`) flag for ignoring the
                           type assumption that polygons are bounded

### Output

`true` if `use_type_assumption` is activated.
Otherwise, `true` iff `P` is bounded.

### Algorithm

If `!use_type_assumption`, we use [`_isbounded_unit_dimensions`](@ref).
"""
function isbounded(P::AbstractHPolygon, use_type_assumption::Bool=true)
    if use_type_assumption
        return true
    end
    if isempty(constraints_list(P))
        return false
    end
    return _isbounded_unit_dimensions(P)
end
