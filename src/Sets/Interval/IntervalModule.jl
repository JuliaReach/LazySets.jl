module IntervalModule

using Reexport

using ..LazySets: AbstractHyperrectangle, LazySet, HalfSpace, Zonotope, UnionSet,
                  _witness_result_empty
using ..LazySets.EmptySetModule: EmptySet
using ..API: eltype, isconvextype, isempty
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Arrays
using ReachabilityBase.Comparison
using ReachabilityBase.Distribution: reseed!
import IntervalArithmetic as IA

@reexport import ..API: an_element, center, constraints_list, diameter, dim, high, ∈,
                        isoperationtype, linear_map, low, permute, rand, rectify, reflect,
                        scale, ρ, σ, translate, vertices_list,
                        difference, intersection, isdisjoint, ⊆, minkowski_difference,
                        minkowski_sum
@reexport import ..LazySets: chebyshev_center_radius, isflat, ngens, plot_recipe,
                             radius_hyperrectangle, split
import Base: convert, -, *, min, max

export Interval

include("Interval.jl")

include("an_element.jl")
include("center.jl")
include("chebyshev_center_radius.jl")
include("constraints_list.jl")
include("convert.jl")
include("diameter.jl")
include("dim.jl")
include("high.jl")
include("in.jl")
include("isflat.jl")
include("linear_map.jl")
include("low.jl")
include("isoperationtype.jl")
include("ngens.jl")
include("permute.jl")
include("radius.jl")
include("radius_hyperrectangle.jl")
include("rand.jl")
include("rectify.jl")
include("reflect.jl")
include("scale.jl")
include("split.jl")
include("support_vector.jl")
include("support_function.jl")
include("translate.jl")
include("vertices_list.jl")

include("difference.jl")
include("intersection.jl")
include("isdisjoint.jl")
include("issubset.jl")
include("minkowski_difference.jl")
include("minkowski_sum.jl")

"""
    -(x::Interval, y::Interval)

Return the difference of two intervals (in the interval-arithmetic sense).

### Input

- `x` -- interval
- `y` -- interval

### Output

The difference of the intervals as a new `Interval` set.
"""
function -(x::Interval, y::Interval)
    return Interval(x.dat - y.dat)
end

"""
```
    *(x::Interval, y::Interval)
```

Return the product of two intervals (in the interval-arithmetic sense).

### Input

- `x` -- interval
- `y` -- interval

### Output

The product of the intervals as a new `Interval` set.
"""
function *(x::Interval, y::Interval)
    return Interval(x.dat * y.dat)
end

"""
    min(x::Interval)

Return the lower component of an interval.

### Input

- `x` -- interval

### Output

The lower (`lo`) component of the interval (a number).
"""
function min(x::Interval)
    return x.dat.lo
end

"""
    max(x::Interval)

Return the higher or upper component of an interval.

### Input

- `x` -- interval

### Output

The higher (`hi`) component of the interval (a number).
"""
function max(x::Interval)
    return x.dat.hi
end

"""
    plot_recipe(x::Interval{N}, [ε]=zero(N)) where {N}

Convert an interval to a pair `(x, y)` of points for plotting.

### Input

- `x` -- interval
- `ε` -- (optional, default: `0`) ignored, used for dispatch

### Output

A pair `(x, y)` of two points that can be plotted.

### Notes

We consider the interval as a line segment with y coordinate equal to zero.
"""
function plot_recipe(x::Interval{N}, ε=zero(N)) where {N}
    return [min(x), max(x)], zeros(N, 2)
end

end  # module
