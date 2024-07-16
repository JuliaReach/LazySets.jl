"""
    constraints_list(x::Interval)

Return the list of constraints of an interval.

### Input

- `x` -- interval

### Output

The list of constraints of the interval represented as two one-dimensional
half-spaces.
"""
function constraints_list(x::Interval)
    require(@__MODULE__, :LazySets; fun_name="constraints_list")

    N = eltype(x)
    constraints = Vector{HalfSpace{N,SingleEntryVector{N}}}(undef, 2)
    e₁ = SingleEntryVector(1, 1, one(N))
    @inbounds constraints[1] = HalfSpace(e₁, max(x))
    @inbounds constraints[2] = HalfSpace(-e₁, -min(x))
    return constraints
end
