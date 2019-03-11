#= concrete implementations of binary intersections between sets =#

export intersection
"""
    intersection(S::AbstractSingleton{N},
                 X::LazySet{N}
                )::Union{Singleton{N}, EmptySet{N}} where {N<:Real}

Return the intersection of a singleton with another set.

### Input

- `S` -- singleton
- `X` -- another set

### Output

If the sets intersect, the result is `S`.
Otherwise, the result is the empty set.
"""
function intersection(S::AbstractSingleton{N},
                      X::LazySet{N}
                     )::Union{Singleton{N}, EmptySet{N}} where {N<:Real}
    return element(S) ∈ X ? S : EmptySet{N}()
end

# symmetric method
function intersection(X::LazySet{N},
                      S::AbstractSingleton{N}
                     )::Union{Singleton{N}, EmptySet{N}} where {N<:Real}
    return intersection(S, X)
end

# disambiguation
function intersection(S1::AbstractSingleton{N},
                      S2::AbstractSingleton{N}
                     )::Union{Singleton{N}, EmptySet{N}} where {N<:Real}
    return element(S1) == element(S2) ? S1 : EmptySet{N}()
end

"""
    intersection(L1::Line{N}, L2::Line{N}
                )::Union{Singleton{N}, Line{N}, EmptySet{N}} where {N<:Real}

Return the intersection of two 2D lines.

### Input

- `L1` -- first line
- `L2` -- second line

### Output

If the lines are identical, the result is the first line.
If the lines are parallel and not identical, the result is the empty set.
Otherwise the result is the only intersection point.

### Examples

The line ``y = -x + 1`` intersected with the line ``y = x``:

```jldoctest
julia> intersection(Line([-1., 1.], 0.), Line([1., 1.], 1.))
Singleton{Float64,Array{Float64,1}}([0.5, 0.5])

julia> intersection(Line([1., 1.], 1.), Line([1., 1.], 1.))
Line{Float64,Array{Float64,1}}([1.0, 1.0], 1.0)
```
"""
function intersection(L1::Line{N}, L2::Line{N}
                     )::Union{Singleton{N}, Line{N}, EmptySet{N}} where {N<:Real}
    b = [L1.b, L2.b]
    a = [transpose(L1.a); transpose(L2.a)]
    try
        # results in LAPACKException or SingularException if parallel
        return Singleton(a \ b)
    catch e
        @assert e isa LAPACKException || e isa SingularException "unexpected " *
            "$(typeof(e)) from LAPACK occurred while intersecting lines:\n$e"
        # lines are parallel
        if an_element(L1) ∈ L2
            # lines are identical
            return L1
        else
            # lines are parallel but not identical
            return EmptySet{N}()
        end
    end
end

"""
    intersection(H1::AbstractHyperrectangle{N},
                 H2::AbstractHyperrectangle{N}
                )::Union{<:Hyperrectangle{N}, EmptySet{N}} where {N<:Real}

Return the intersection of two hyperrectangles.

### Input

- `H1` -- first hyperrectangle
- `H2` -- second hyperrectangle

### Output

If the hyperrectangles do not intersect, the result is the empty set.
Otherwise the result is the hyperrectangle that describes the intersection.

### Algorithm

In each isolated direction `i` we compute the rightmost left border and the
leftmost right border of the hyperrectangles.
If these borders contradict, then the intersection is empty.
Otherwise the result uses these borders in each dimension.
"""
function intersection(H1::AbstractHyperrectangle{N},
                      H2::AbstractHyperrectangle{N}
                     )::Union{Hyperrectangle{N}, EmptySet{N}} where {N<:Real}
    n = dim(H1)
    c1 = center(H1)
    c2 = center(H2)
    r1 = radius_hyperrectangle(H1)
    r2 = radius_hyperrectangle(H2)
    high = Vector{N}(undef, n)
    low = Vector{N}(undef, n)
    for i in 1:n
        high1 = c1[i] + r1[i]
        low1 = c1[i] - r1[i]
        high2 = c2[i] + r2[i]
        low2 = c2[i] - r2[i]
        high[i] = min(high1, high2)
        low[i] = max(low1, low2)
        if high[i] < low[i]
            return EmptySet{N}()
        end
    end
    return Hyperrectangle(high=high, low=low)
end

# disambiguation
function intersection(S::AbstractSingleton{N},
                      H::AbstractHyperrectangle{N}
                     )::Union{Singleton{N}, EmptySet{N}} where {N<:Real}
    return invoke(intersection, Tuple{typeof(S), LazySet{N}}, S, H)
end
function intersection(H::AbstractHyperrectangle{N},
                      S::AbstractSingleton{N}
                     )::Union{Singleton{N}, EmptySet{N}} where {N<:Real}
    return invoke(intersection, Tuple{typeof(S), LazySet{N}}, S, H)
end

"""
    intersection(x::Interval{N},
                 y::Interval{N}
                 )::Union{Interval{N}, EmptySet{N}} where {N<:Real}

Return the intersection of two intervals.

### Input

- `x` -- first interval
- `y` -- second interval

### Output

If the intervals do not intersect, the result is the empty set.
Otherwise the result is the interval that describes the intersection.
"""
function intersection(x::Interval{N},
                      y::Interval{N}
                     )::Union{Interval{N}, EmptySet{N}} where {N<:Real}
    if min(y) > max(x) || min(x) > max(y)
        return EmptySet{N}()
    else
        return Interval(max(min(x), min(y)), min(max(x), max(y)))
    end
end

"""
    intersection(P1::AbstractHPolygon{N},
                 P2::AbstractHPolygon{N}
                )::Union{HPolygon{N}, EmptySet{N}} where {N<:Real}

Return the intersection of two polygons in constraint representation.

### Input

- `P1`    -- first polygon
- `P2`    -- second polygon
- `prune` -- (optional, default: `true`) flag for removing redundant constraints

### Output

If the polygons do not intersect, the result is the empty set.
Otherwise the result is the polygon that describes the intersection.

### Algorithm

We just combine the constraints of both polygons.
To obtain a linear-time algorithm, we interleave the constraints.
If there are two constraints with the same normal vector, we choose the tighter
one.

Redundancy of constraints is checked with
[`remove_redundant_constraints!(::AbstractHPolygon)`](@ref).
"""
function intersection(P1::AbstractHPolygon{N},
                      P2::AbstractHPolygon{N},
                      prune::Bool=true
                     )::Union{HPolygon{N}, EmptySet{N}} where {N<:Real}
    # all constraints of one polygon are processed; now add the other polygon's
    # constraints
    @inline function add_remaining_constraints!(c, i, c1, i1, duplicates)
        c[i+1:length(c)-duplicates] = c1[i1:length(c1)]
        return true
    end

    # choose the constraint of c1; the directions are different
    @inline function choose_first_diff_dir!(c, i, i1, i2, c1, c2, duplicates)
        c[i] = c1[i1]
        if i1 == length(c1)
            return add_remaining_constraints!(c, i, c2, i2, duplicates)
        end
        return false
    end

    # choose the constraint of c1; the directions are equivalent (i.e., linearly
    # dependent)
    @inline function choose_first_same_dir!(c, i, i1, i2, c1, c2, duplicates)
        c[i] = c1[i1]
        if i1 == length(c1)
            if i2 < length(c2)
                return add_remaining_constraints!(c, i, c2, i2+1, duplicates)
            end
            return true
        elseif i2 == length(c2)
            return add_remaining_constraints!(c, i, c1, i1+1, duplicates)
        end
        return false
    end

    c1 = constraints_list(P1)
    c2 = constraints_list(P2)
    if length(c1) == 0
        return P2
    elseif length(c2) == 0
        return P1
    end
    c = Vector{LinearConstraint{N}}(undef, length(c1) + length(c2))
    i1 = 1
    i2 = 1
    duplicates = 0
    for i in 1:length(c)
        if c1[i1].a <= c2[i2].a
            if c2[i2].a <= c1[i1].a
                duplicates += 1
                # constraints have the same normal vector: take the tighter one
                if is_tighter_same_dir_2D(c1[i1], c2[i2])
                    # first constraint is tighter
                    if choose_first_same_dir!(c, i, i1, i2, c1, c2, duplicates)
                        break
                    end
                else
                    # second constraint is tighter
                    if choose_first_same_dir!(c, i, i2, i1, c2, c1, duplicates)
                        break
                    end
                end
                i1 += 1
                i2 += 1
            else
                # first constraint comes first
                if choose_first_diff_dir!(c, i, i1, i2, c1, c2, duplicates)
                    break
                end
                i1 += 1
            end
        else
            # second constraint comes first
            if choose_first_diff_dir!(c, i, i2, i1, c2, c1, duplicates)
                break
            end
            i2 += 1
        end
    end
    if duplicates > 0
        deleteat!(c, length(c)-duplicates+1:length(c))
    end

    P = HPolygon(c, sort_constraints=false)
    if prune
        remove_redundant_constraints!(P)
        if isempty(P)
            return EmptySet{N}()
        end
    end
    return P
end

"""
    intersection(P1::AbstractPolyhedron{N},
                 P2::AbstractPolyhedron{N};
                 backend=nothing,
                 use_polyhedra_interface=false) where {N<:Real}

Compute the intersection of two polyhedra.

### Input

- `P1`        -- polyhedron
- `P2`        -- polyhedron
- `backend`   -- (optional, default: `nothing`) the LP solver or the backend for
                 polyhedral computations; its value is set internally, see the
                 Notes below for details
- `use_polyhedra_interface` -- (optional, default: `false`) if `true`, use the
                 `Polyhedra` interface for the removal of constraints

### Output

An `HPolyhedron` resulting from the intersection of `P1` and `P2`, with the
redundant constraints removed, or an empty set if the intersection is empty.
If one of the arguments is a polytope, the result is an `HPolytope` instead.

### Notes

The default value of the backend is set internally and depends on whether the
Polyhedra backend is used or not. The default backends are `GLPKSolverLP()`
and `default_polyhedra_backend(P1, N)`, respectively.

Note that if `use_polyhedra_interface` is set to `true`, there is no guarantee
that the removal of constraints keep the set empty (see #1038 and
Polyhedra#146), so it is better to check for emptiness of intersection before
using this function in that case.

### Algorithm

This implementation unifies the constraints of the two sets obtained from the
`constraints_list` method.
"""
function intersection(P1::AbstractPolyhedron{N},
                      P2::AbstractPolyhedron{N};
                      backend=nothing,
                      use_polyhedra_interface=false) where {N<:Real}
    HPOLY = (P1 isa AbstractPolytope || P2 isa AbstractPolytope) ?
        HPolytope{N} : HPolyhedron{N}

    # concatenate the linear constraints
    Q = HPOLY([constraints_list(P1); constraints_list(P2)])

    # remove redundant constraints
    if use_polyhedra_interface
        if backend == nothing
            backend = default_polyhedra_backend(P1, N)
        end
        # convert to an hrep, remove the redundancies and convert back to HPOLY
        ph = polyhedron(Q; backend=backend)
        removehredundancy!(ph)
        return convert(HPOLY, ph)
    else
        if backend == nothing
            backend = GLPKSolverLP()
        end
        # here, detection of empty intersection may be reported as an infeasible LP
        if remove_redundant_constraints!(Q, backend=backend)
            return Q
        else
            return EmptySet{N}()
        end
    end
end

# disambiguation
function intersection(S::AbstractSingleton{N},
                      P::AbstractPolyhedron{N}
                     )::Union{Singleton{N}, EmptySet{N}} where {N<:Real}
    return invoke(intersection, Tuple{typeof(S), LazySet{N}}, S, P)
end
function intersection(P::AbstractPolyhedron{N},
                      S::AbstractSingleton{N}
                     )::Union{Singleton{N}, EmptySet{N}} where {N<:Real}
    return invoke(intersection, Tuple{typeof(S), LazySet{N}}, S, P)
end

"""
    intersection(P1::Union{VPolytope{N}, VPolygon{N}},
                 P2::Union{VPolytope{N}, VPolygon{N}};
                 [backend]=default_polyhedra_backend(P1, N),
                 [prunefunc]=removevredundancy!) where {N<:Real}

Compute the intersection of two polytopes in vertex representation.

### Input

- `P1`        -- polytope in vertex representation
- `P2`        -- polytope in vertex representation
- `backend`   -- (optional, default: `default_polyhedra_backend(P1, N)`) the
                 backend for polyhedral computations
- `prunefunc` -- (optional, default: `removevredundancy!`) function to prune
                 the vertices of the result

### Output

A `VPolygon` if both arguments are `VPolygon`s, and a `VPolytope` otherwise.
"""
function intersection(P1::Union{VPolytope{N}, VPolygon{N}},
                      P2::Union{VPolytope{N}, VPolygon{N}};
                      backend=default_polyhedra_backend(P1, N),
                      prunefunc=removevredundancy!) where {N<:Real}
    Q1 = polyhedron(convert(VPolytope, P1); backend=backend)
    Q2 = polyhedron(convert(VPolytope, P2); backend=backend)
    Pint = Polyhedra.intersect(Q1, Q2)
    prunefunc(Pint)
    res = VPolytope(Pint)
    if P1 isa VPolygon && P2 isa VPolygon
        return convert(VPolygon, res)
    end
    return res
end

"""
    intersection(cup::UnionSet{N}, X::LazySet{N}) where {N<:Real}

Return the intersection of a union of two convex sets and another convex set.

### Input

- `cup` -- union of two convex sets
- `X`   -- convex set

### Output

The union of the pairwise intersections, expressed as a `UnionSet`.
If one of those sets is empty, only the other set is returned.
"""
function intersection(cup::UnionSet{N}, X::LazySet{N}) where {N<:Real}
    return intersection(cup.X, X) ∪ intersection(cup.Y, X)
end

# symmetric method
function intersection(X::LazySet{N}, cup::UnionSet{N}) where {N<:Real}
    return intersection(cup, X)
end

# disambiguation
function intersection(cup::UnionSet{N}, S::AbstractSingleton{N}) where {N<:Real}
    return element(S) ∈ cup ? S : EmptySet{N}()
end
function intersection(S::AbstractSingleton{N}, cup::UnionSet{N}) where {N<:Real}
    return invoke(intersection, Tuple{UnionSet{N}, typeof(S)}, cup, S)
end

"""
    intersection(cup::UnionSetArray{N}, X::LazySet{N}) where {N<:Real}

Return the intersection of a union of a finite number of convex sets and another
convex set.

### Input

- `cup` -- union of a finite number of convex sets
- `X`   -- convex set

### Output

The union of the pairwise intersections, expressed as a `UnionSetArray`.
"""
function intersection(cup::UnionSetArray{N}, X::LazySet{N}) where {N<:Real}
    return UnionSetArray([intersection(Y, X) for Y in array(cup)])
end

# symmetric method
function intersection(X::LazySet{N}, cup::UnionSetArray{N}) where {N<:Real}
    return intersection(cup, X)
end

# disambiguation
function intersection(cup::UnionSetArray{N},
                      S::AbstractSingleton{N}) where {N<:Real}
    return element(S) ∈ cup ? S : EmptySet{N}()
end
function intersection(S::AbstractSingleton{N},
                      cup::UnionSetArray{N}) where {N<:Real}
    return invoke(intersection, Tuple{UnionSetArray{N}, typeof(S)}, cup, S)
end

"""
    intersection(L::LinearMap{N}, S::LazySet{N}) where {N}

Return the intersection of a lazy linear map and a convex set.

### Input

- `L` -- linear map
- `S` -- convex set
  
### Output

The polytope obtained by the intersection of `l.M * L.X` and `S`.
"""
function intersection(L::LinearMap{N}, S::LazySet{N}) where {N}
    return intersection(linear_map(L.M, L.X), S)
end

# symmetric method
function intersection(S::LazySet{N}, L::LinearMap{N}) where {N}
    return intersection(L, S)
end

# disambiguation
function intersection(L1::LinearMap{N}, L2::LinearMap{N}) where {N}
    return intersection(linear_map(L1.M, L1.X), linear_map(L2.M, L2.X))
end

"""
function intersection(X::CartesianProductArray{N},
                      Y::AbstractPolyhedron{N}) where {N}

Return the intersection of the cartesian product of a finite number of convex sets and a polyhedron.

### Input

 - `X` -- cartesian product of a finite number of convex sets
 - `Y` -- polyhedron

### Output

The concrete intersection between `X` and `Y`

### Algorithm

This function takes into account that if the polyhedron is unbounded, the intersection
is only needed to be taken in the elements of the cartesian product array (subsets of
variables, or "blocks") which are constrained; those which are not constrained do not need to be intersected.
"""

function intersection_blocks(X::CartesianProductArray{N},
                      Y::AbstractPolyhedron{N}, blocks::Dict{Int,Int}) where {N}
    # preallocate the resulting CartesianProductArray with size hint
    result = CartesianProductArray(length(X.array), N)

    for bi in 1:length(X.array)
        if haskey(blocks,bi)
            # otherwise, make the intersection with the projection of the halfspace
            vars = variable_indices(X, bi, blocks[bi])
            push!(result.array, intersection(X.array[bi], Approximations.project(Y,vars, LinearMap)))
        else
            # if this block is not constrained, just push the set
            push!(result.array, X.array[bi])
        end
    end
    return result
end

function intersection_combine(X::CartesianProductArray{N},
                      Y::AbstractPolyhedron{N}, blocks::Dict{Int,Int}) where {N}
    result = CartesianProductArray(length(X.array), N)

    low_set = CartesianProductArray(length(blocks), N)
    vars = Vector{Int}()
    block_structure = Vector{Int}()
    blocks = sort(blocks)

    for bi in keys(blocks)
        push!(low_set.array, X.array[bi])
        append!(vars, variable_indices(X, bi, blocks[bi]))
        push!(block_structure, dim(X.array[bi]))
    end

    approx_low_set = HPolytope(constraints_list(low_set));

    low_intersection = intersection(approx_low_set, Approximations.project(Y,vars, LinearMap))

    if isempty(low_intersection)
        return EmptySet()
    end
    decomposed_low_set = Approximations.decompose(low_intersection, set_type=LinearMap, blocks=block_structure)


    index = 1
    for bi in 1:length(X.array)
        if haskey(blocks,bi)
            push!(result.array, decomposed_low_set.array[index])
            index += 1
        else
            push!(result.array, X.array[bi])
        end
    end
    return result
end

function intersection(X::CartesianProductArray{N},
                      Y::AbstractPolyhedron{N},
                      blocks::Dict{Int,Int},
                      is_combine::Bool) where {N}
    if is_combine
        return intersection_combine(X, Y, blocks)
    else
        return intersection_blocks(X, Y, blocks)
    end
end

function intersection(X::CartesianProductArray{N},
                      Y::AbstractPolyhedron{N}, is_combine::Bool=true) where {N}

    if isbounded(Y)
        # no free variables
        blocks = block_indices(X)
    else
        constrained_vars = constrained_dimensions(Y)
        blocks = block_indices(X, constrained_vars)
    end
    return intersection(X, Y, blocks, is_combine)
end


# symmetric method
function intersection(Y::AbstractPolyhedron{N}, X::CartesianProductArray{N}, is_combine::Bool=true) where {N}
    intersection(X, Y, is_combine)
end


function intersection(X::CartesianProductArray{N},
                      Y::AbstractPolyhedron{N}, vars::AbstractVector) where {N}
    # preallocate the resulting CartesianProductArray with size hint
    result = CartesianProductArray(length(X.array), N)

    blocks = block_indices(X, vars)

    for bi in 1:length(X.array)
        if haskey(blocks,bi)
            # otherwise, make the intersection with the projection of the halfspace
            vars = variable_indices(X, bi, blocks[bi])
            push!(result.array,intersection(X.array[bi], Approximations.project(Y,vars, LinearMap)))
        end
    end
    return result
end

"""
    intersection(U::Universe{N}, X::LazySet{N}) where {N<:Real}

Return the intersection of a universe and a convex set.

### Input

- `U` -- universe
- `X` -- convex set

### Output

The set `X`.
"""
function intersection(U::Universe{N}, X::LazySet{N}) where {N<:Real}
    return X
end

# symmetric method
function intersection(X::LazySet{N}, U::Universe{N}) where {N<:Real}
    return X
end

# disambiguation
function intersection(U::Universe{N}, ::Universe{N}) where {N<:Real}
    return U
end
function intersection(U::Universe{N}, P::AbstractPolyhedron{N}) where {N<:Real}
    return P
end
function intersection(P::AbstractPolyhedron{N}, U::Universe{N}) where {N<:Real}
    return P
end
function intersection(U::Universe{N}, S::AbstractSingleton{N}) where {N<:Real}
    return S
end
function intersection(S::AbstractSingleton{N}, U::Universe{N}) where {N<:Real}
    return S
end

"""
    intersection(P::AbstractPolyhedron{N}, rm::ResetMap{N}) where {N<:Real}

Return the intersection of a polyhedron and a polyhedral reset map.

### Input

- `P`  -- polyhedron
- `rm` -- polyhedral reset map

### Output

A polyhedron.

### Notes

We assume that `rm` is polyhedral, i.e., has a `constraints_list` method
defined.
"""
function intersection(P::AbstractPolyhedron{N}, rm::ResetMap{N}) where {N<:Real}
    return intersection(P, HPolyhedron(constraints_list(rm)))
end

# symmetric method
function intersection(rm::ResetMap{N}, P::AbstractPolyhedron{N}) where {N<:Real}
    return intersection(P, rm)
end

# more efficient version for polytopic
function intersection(P::AbstractPolyhedron{N},
                      rm::ResetMap{N, <:AbstractPolytope}) where {N<:Real}
    return intersection(P, HPolytope(constraints_list(rm)))
end

# symmetric method
function intersection(rm::ResetMap{N, <:AbstractPolytope},
                      P::AbstractPolyhedron{N}) where {N<:Real}
    return intersection(P, rm)
end
