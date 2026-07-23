import LazySets

using LazySets: @validate
using LazySets.EmptySetModule: EmptySet
using LazySets.HalfSpaceModule: HalfSpace
using LazySets.IntervalModule: Interval, _center
using LazySets.ZonotopeModule: Zonotope
using LazySets: UnionSet
using ReachabilityBase.Arrays: SingleEntryVector
using ReachabilityBase.Comparison: isapproxzero
import LazySets.IntervalModule: _constraints_list_Vector, _linear_map_zonotope
import LazySets.API: complement, constraints_list, difference, intersection,
                     minkowski_difference

function complement(X::Interval)
    N = eltype(X)
    L = HalfSpace(SingleEntryVector(1, 1, one(N)), min(X))
    H = HalfSpace(SingleEntryVector(1, 1, -one(N)), -max(X))
    return UnionSet(L, H)
end

function constraints_list(X::Interval)
    N = eltype(X)
    constraints = Vector{HalfSpace{N,SingleEntryVector{N}}}(undef, 2)
    e₁ = SingleEntryVector(1, 1, one(N))
    @inbounds constraints[1] = HalfSpace(e₁, max(X))
    @inbounds constraints[2] = HalfSpace(-e₁, -min(X))
    return constraints
end

function _constraints_list_Vector(X::Interval)
    N = eltype(X)
    constraints = Vector{HalfSpace{N,Vector{N}}}(undef, 2)
    @inbounds constraints[1] = HalfSpace([one(N)], max(X))
    @inbounds constraints[2] = HalfSpace([-one(N)], -min(X))
    return constraints
end

"""
# Extended help

    difference(X::Interval, Y::Interval)

### Output

Depending on the position of the intervals, the output is one of the following:

- An `EmptySet`.
- An `Interval`.
- A `UnionSet` of two `Interval` sets.

### Algorithm

Let ``X = [a, b]`` and ``Y = [c, d]`` be intervals. Their set difference is
``X ∖ Y = \\{x: x ∈ X \\text{ and } x ∉ Y \\}`` and, depending on their
position, three different results may occur:

- If ``X`` and ``Y`` do not overlap, i.e., if their intersection is empty, then
  the set difference is just ``X``.
- Otherwise, let ``Z = X ∩ Y ≠ ∅``, then ``Z`` splits ``X`` into either one or
  two intervals. The latter case happens when the bounds of ``Y`` are strictly
  contained in ``X``.

To check for strict inclusion, we assume that the inclusion is strict and then
check whether the resulting intervals that cover ``X`` (one to its left and one
to its right, let them be `L` and `R`), obtained by intersection with ``Y``, are
flat. Three cases may arise:

- If both `L` and `R` are flat then ``X = Y`` and the result is the empty set.
- If only `L` is flat, then the result is `R`, the remaining interval not
  covered by ``Y``. Similarly, if only `R` is flat, then the result is `L`.
- Finally, if none of the intervals is flat, then ``Y`` is strictly contained
  in ``X`` and the set union of `L` and `R` is returned.

### Examples

```@meta
DocTestSetup = quote
    using LazySets
end
```

```jldoctest
julia> X = Interval(0, 2); Y = Interval(1, 4); Z = Interval(2, 3);

julia> difference(X, X)
∅(1)

julia> difference(X, Y)
Interval{Float64}([0.0, 1.0])

julia> difference(Y, Z)
UnionSet{Float64, Interval{Float64}, Interval{Float64}}(Interval{Float64}([1.0, 2.0]), Interval{Float64}([3.0, 4.0]))
```
"""
@validate function difference(X::Interval, Y::Interval)
    l, h = _intersection_interval_bounds(X, Y)
    if l > h
        return X
    else
        flat_left = isapproxzero(l - min(X))
        flat_right = isapproxzero(max(X) - h)

        if flat_left && flat_right
            N = promote_type(eltype(X), eltype(Y))
            return EmptySet{N}(1)
        elseif flat_left && !flat_right
            return Interval(h, max(X))
        elseif !flat_left && flat_right
            return Interval(min(X), l)
        else
            return UnionSet(Interval(min(X), l), Interval(h, max(X)))
        end
    end
end

@validate function intersection(X::Interval, Y::Interval)
    l, h = _intersection_interval_bounds(X, Y)
    if l > h
        N = promote_type(eltype(X), eltype(Y))
        return EmptySet{N}(1)
    else
        return Interval(l, h)
    end
end

function _intersection_interval_bounds(X::Interval, Y::Interval)
    l = max(min(X), min(Y))
    h = min(max(X), max(Y))
    return (l, h)
end

function _linear_map_zonotope(M::AbstractMatrix, X::Interval)
    nout = size(M, 1)
    cx = _center(X)
    gx = cx - min(X)
    N = promote_type(eltype(M), eltype(X))
    c = Vector{N}(undef, nout)
    gen = Matrix{N}(undef, nout, 1)
    @inbounds for i in 1:nout
        c[i] = M[i, 1] * cx
        gen[i] = M[i, 1] * gx
    end
    return Zonotope(c, gen)
end

"""
# Extended help

    minkowski_difference(I1::Interval, I2::Interval)

### Output

An `Interval`, or an `EmptySet` if the difference is empty.
"""
@validate function minkowski_difference(I1::Interval, I2::Interval)
    l = min(I1) - min(I2)
    h = max(I1) - max(I2)
    if h < l
        N = promote_type(eltype(I1), eltype(I2))
        return EmptySet{N}(1)
    end
    return Interval(l, h)
end
