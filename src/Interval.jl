@require IntervalArithmetic begin

import IntervalArithmetic.AbstractInterval

"""
    IA = IntervalArithmetic

Name alias for the interval arithmetic library.
"""
const IA = IntervalArithmetic

export Interval, IA,
       dim, σ, center, +, -, *, low, high, vertices_list

"""
    Interval{N, IN <: IA.AbstractInterval{N}} <: LazySet{N}

Type representing an interval on the real line.

### Fields

- `dat` -- data container for the given interval

### Examples

```jldoctest interval_constructor
julia> using LazySets, IntervalArithmetic

julia> x = LazySets.Interval(IA.Interval(0.0, 1.0))
LazySets.Interval{Float64,IntervalArithmetic.Interval{Float64}}([0, 1])

julia> dim(x)
1

julia> center(x)
0.5
```
"""
struct Interval{N, IN <: AbstractInterval{N}} <: AbstractPointSymmetricPolytope{N}
   dat::IN
end
# type-less convenience constructor
Interval{N, IN <: AbstractInterval{N}}(interval::IN) = Interval{N, IN}(interval)
# TODO: adapt show method

"""
    dim(x::Interval)::Int

Return the ambient dimension of an interval.

### Input

- `x` -- interval

### Output

The integer 1.
"""
dim(x::Interval)::Int = 1

"""
    σ(d::AbstractVector{N}, x::Interval{N, IN})::AbstractVector{N} where {N, IN <: IA.AbstractInterval{N}}

Return the support vector of an ellipsoid in a given direction.

### Input

- `d` -- direction
- `x` -- interval

### Output

Support vector in the given direction.
"""
function σ(d::AbstractVector{N}, x::Interval{N, IN})::AbstractVector{N} where {N, IN <: IA.AbstractInterval{N}}
    return d[1] > 0 ? [x.dat.hi] : [x.dat.lo]
end

import LazySets.center

"""
    center(x::Interval)

Return the center of ther interval.

### Input

- `x` -- interval

### Output

The center, or midpoint, of ``x``.
"""
center(x::Interval) = IA.mid(x.dat)

import Base:+, -, *

+(x::Interval, y::Interval) = Interval(x.dat + y.dat)
-(x::Interval, y::Interval) = Interval(x.dat - y.dat)
*(x::Interval, y::Interval) = Interval(x.dat * y.dat)

"""
    low(x::Interval)

Return the lower component of an interval.

### Input

- `x` -- interval

### Output

The lower (`lo`) component of the interval.
"""
low(x::Interval) = x.dat.lo

"""
    high(x::Interval)

Return the higher or upper component of an interval.

### Input

- `x` -- interval

### Output

The higher (`hi`) component of the interval.
"""
high(x::Interval) = x.dat.hi


"""
    vertices_list(x::Interval)

Return the list of vertices of this interval.

### Input

- `x` -- interval

### Output

The list of vertices of the interval represented as two one-dimensional vectors.
"""
vertices_list(x::Interval) = [[low(x::Interval)], [high(x::Interval)]]

end
