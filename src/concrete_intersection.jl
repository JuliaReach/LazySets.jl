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

- `P1` -- first polygon
- `P2` -- second polygon

### Output

If the polygons do not intersect, the result is the empty set.
Otherwise the result is the polygon that describes the intersection.

### Algorithm

We just combine the constraints of both polygons.
To obtain a linear-time algorithm, we interleave the constraints.
If there are two constraints with the same normal vector, we choose the tighter
one.
"""
function intersection(P1::AbstractHPolygon{N},
                      P2::AbstractHPolygon{N}
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

    # check if the first of two constraint with equivalent direction is tighter
    @inline function is_first_constraint_tighter(lc1::LinearConstraint{N},
                                                 lc2::LinearConstraint{N}
                                                ) where {N<:Real}
        if lc1.a[1] == zero(N)
            @assert lc2.a[1] == zero(N)
            return lc1.b <= lc1.a[2]/lc2.a[2] * lc2.b
        end
        return lc1.b <= lc1.a[1]/lc2.a[1] * lc2.b
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
                if is_first_constraint_tighter(c1[i1], c2[i2])
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

    # TODO: remove redundant constraints (#582) and return an EmptySet if empty
    return P
end

"""
    intersection(P::HPoly{N},
                 hs::HalfSpace{N};
                 backend=default_polyhedra_backend(P, N),
                 prunefunc=removehredundancy!) where {N<:Real}

Compute the intersection of a polytope in H-representation and a half-space.

### Input

- `P`         -- polytope
- `hs`        -- half-space
- `backend`   -- (optional, default: `nothing`) the LP solver or the the
                 polyhedral computations backend; its value is set internally,
                 see see below in the Notes for details
- `use_polyhedra_interface` -- (optional, default: `false`) if `true`, use the
                 `Polyhedra` interface for the removal of constraints
### Output

The same polytope in H-representation with just one more constraint.
"""
function intersection(P::HPoly{N},
                      hs::HalfSpace{N};
                      backend=nothing,
                      use_polyhedra_interface=false) where {N<:Real}
    return intersection(P, HPolyhedron([hs]), backend=backend, use_polyhedra_interface=use_polyhedra_interface)
end

# symmetric method
function intersection(hs::HalfSpace{N},
                      P::HPoly{N};
                      backend=nothing,
                      use_polyhedra_interface=false) where {N<:Real}
    return intersection(P, hs, backend=backend, use_polyhedra_interface=use_polyhedra_interface)
end

"""
    intersection(P1::HPoly{N},
                 P2::HPoly{N};
                 backend=nothing,
                 use_polyhedra_interface=false) where {N<:Real}

Compute the intersection of two polyhedra in H-representation.

### Input

- `P1`        -- polyhedron
- `P2`        -- polyhedron
- `backend`   -- (optional, default: `nothing`) the LP solver or the the
                 polyhedral computations backend; its value is set internally,
                 see see below in the Notes for details
- `use_polyhedra_interface` -- (optional, default: `false`) if `true`, use the
                 `Polyhedra` interface for the removal of constraints

### Output

A polyhedron resulting from the intersection of `P1` and `P2`, with the redundant
constraints removed, or an empty set if the intersection is empty.

### Notes

The default value of the backend is set internally and depends on whether the
Polyhedra backend is used or not. The default backends are `GLPKSolverLP()`
and `default_polyhedra_backend(P1, N)` respectively.

Note that if `use_polyhedra_interface` is set to `true`, there is no guarantee
that the removal of constraints keep the set empty (see #1038 and Polyhedra#146),
so it is better to check for emptiness of intersection before using this function
in that case.

The method implemented in this function can be used for any pair of sets that can
handle the `constraints_list` option.
"""
function intersection(P1::HPoly{N},
                      P2::HPoly{N};
                      backend=nothing,
                      use_polyhedra_interface=false) where {N<:Real}
    if typeof(P1) == typeof(P2)
        HPOLY = typeof(P1)
    else
        # one of them must be a polytope, so the intersection will be bounded
        HPOLY = HPolytope{N}
    end

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

"""
    intersection(P1::HPoly{N},
                 P2::VPolytope{N};
                 backend=default_polyhedra_backend(P1, N),
                 prunefunc=removehredundancy!) where {N<:Real}

Compute the intersection of two polytopes in either H-representation or
V-representation.

### Input

- `P1`        -- polytope
- `P2`        -- polytope
- `backend`   -- (optional, default: `default_polyhedra_backend(P1, N)`) the
                 polyhedral computations backend, see
                 [Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/latest/installation.html#Getting-Libraries-1)
                 for further information
- `prunefunc` -- (optional, default: `removehredundancy!`) function to
                 post-process the output of `intersect`

### Output

The polytope obtained by the intersection of `P1` and `P2`.
"""
function intersection(P1::HPoly{N},
                      P2::VPolytope{N};
                      backend=default_polyhedra_backend(P1, N),
                      prunefunc=removehredundancy!) where {N<:Real}

    Q1 = polyhedron(P1; backend=backend)
    Q2 = polyhedron(P2; backend=backend)
    Pint = Polyhedra.intersect(Q1, Q2)
    prunefunc(Pint)
    return convert(typeof(P1), Pint)
end

# symmetric one
function intersection(P1::VPolytope{N},
                      P2::HPoly{N};
                      backend=default_polyhedra_backend(P1, N),
                      prunefunc=removehredundancy!) where {N<:Real}
    return intersection(P2, P1; backend=backend, prunefunc=prunefunc)
end

# missing case: VRep with VRep
function intersection(P1::VPolytope{N},
                      P2::VPolytope{N};
                      backend=default_polyhedra_backend(P1, N),
                      prunefunc=removehredundancy!) where {N<:Real}
    Q1 = polyhedron(P1; backend=backend)
    Q2 = polyhedron(P2; backend=backend)
    Pint = Polyhedra.intersect(Q1, Q2)
    prunefunc(Pint)
    return VPolytope(Pint)
end

"""
    intersection(P1::HPoly{N},
                 P2::AbstractPolytope{N};
                 backend=default_polyhedra_backend(P1, N),
                 prunefunc=removehredundancy!) where {N<:Real}

Compute the intersection of a polyhedron and a polytope.

### Input

- `P1`        -- polyhedron
- `P2`        -- polytope
- `backend`   -- (optional, default: `default_polyhedra_backend(P1, N)`) the
                 polyhedral computations backend, see
                 [Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/latest/installation.html#Getting-Libraries-1)
                 for further information
- `prunefunc` -- (optional, default: `removehredundancy!`) function to
                 post-process the output of `intersect`

### Output

The polytope in H-representation obtained by the intersection of `P1` and `P2`.
"""
function intersection(P1::HPoly{N},
                      P2::AbstractPolytope{N};
                      backend=default_polyhedra_backend(P1, N),
                      prunefunc=removehredundancy!) where {N<:Real}
    return intersection(P1, HPolytope(constraints_list(P2)))
end

# symmetric function
function intersection(P1::AbstractPolytope{N},
                      P2::HPoly{N};
                      backend=default_polyhedra_backend(P1, N),
                      prunefunc=removehredundancy!) where {N<:Real}
    return intersection(P2, P1)
end

# disambiguation
function intersection(S::AbstractSingleton{N},
                      P::HPoly{N}
                     )::Union{Singleton{N}, EmptySet{N}} where {N<:Real}
    return invoke(intersection, Tuple{typeof(S), LazySet{N}}, S, P)
end
function intersection(P::HPoly{N},
                      S::AbstractSingleton{N}
                     )::Union{Singleton{N}, EmptySet{N}} where {N<:Real}
    return invoke(intersection, Tuple{typeof(S), LazySet{N}}, S, P)
end


"""
    intersection(P1::S1, P2::S2) where {S1<:AbstractPolytope{N},
                                        S2<:AbstractPolytope{N}} where {N<:Real}

Compute the intersection of two polytopic sets.

### Input

- `P1` -- polytope
- `P2` -- another polytope

### Output

The polytope obtained by the intersection of `P1` and `P2`.
Usually the V-representation is used.

### Notes

This fallback implementation requires `Polyhedra` to evaluate the concrete
intersection.
Inputs that are not of type `HPolytope` or `VPolytope` are converted to an
`HPolytope` through the `constraints_list` function.
"""
function intersection(P1::S1, P2::S2) where {S1<:AbstractPolytope{N},
                                             S2<:AbstractPolytope{N}} where
                                            {N<:Real}
    function get_polytope(P::Union{HPolytope, VPolytope})
        return P
    end
    function get_polytope(P::AbstractPolytope)
        return HPolytope(constraints_list(P))
    end
    return intersection(get_polytope(P1), get_polytope(P2))
end

# disambiguation
function intersection(S::AbstractSingleton{N},
                      P::AbstractPolytope{N}
                     )::Union{Singleton{N}, EmptySet{N}} where {N<:Real}
    return invoke(intersection, Tuple{typeof(S), LazySet{N}}, S, P)
end

function intersection(P::AbstractPolytope{N},
                      S::AbstractSingleton{N}
                     )::Union{Singleton{N}, EmptySet{N}} where {N<:Real}
    return invoke(intersection, Tuple{typeof(S), LazySet{N}}, S, P)
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
intersection(ca1::CartesianProductArray{N, S}, L2::CartesianProductArray{N, S}, nonzero_blocks::Vector{Int})
              ) where {N<:Real, S<:LazySet{N}}

Return the intersection of two CartesianProductArray's only for necessary blocks.

### Input

 - `ca1` -- Cartesian Product Array of convex sets
 - `ca2` -- Cartesian Product Array of LazySets (it might by unbounded sets)
 - `nonzero_blocks` -- Required blocks to intersect

### Output

The Cartesian Product Array obtained by the intersection nonzero blocks of `ca1` and `ca2`.
"""
function intersection(ca1::CartesianProductArray{N, S}, L2::CartesianProductArray{N, S}, nonzero_blocks::Vector{Int})
              ) where {N<:Real, S<:LazySet{N}}
    for i in nonzero_blocks
        ca1.array[i] = intersection(ca1.array[i], ca2.array[i])
    end
    return ca1
end
