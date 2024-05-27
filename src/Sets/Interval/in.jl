"""
    ∈(v::AbstractVector, x::Interval))

Check whether a given point is contained in an interval.

### Input

- `v` -- point/vector
- `x` -- interval

### Output

`true` iff ``v ∈ x``.
"""
function ∈(v::AbstractVector, x::Interval)
    @assert length(v) == 1 "a $(length(v))-dimensional vector is " *
                           "incompatible with an interval"
    return @inbounds v[1] ∈ x.dat
end

"""
    ∈(v::Number, x::Interval)

Check whether a number is contained in an interval.

### Input

- `v` -- scalar
- `x` -- interval

### Output

`true` iff `x` contains the singleton `[v]`.
"""
function ∈(v::Number, x::Interval)
    return v ∈ x.dat
end
