@require IntervalArithmetic begin

import IntervalArithmetic.AbstractInterval

"""
    IA = IntervalArithmetic

Name alias for the interval arithmetic library.
"""
const IA = IntervalArithmetic

export Interval, IA

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

end
