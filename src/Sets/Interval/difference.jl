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

```jldoctest
julia> X = Interval(0, 2); Y = Interval(1, 4); Z = Interval(2, 3);

julia> difference(X, X)
∅(1)

julia> difference(X, Y)
Interval{Float64}([0, 1])

julia> difference(Y, Z)
UnionSet{Float64, Interval{Float64}, Interval{Float64}}(Interval{Float64}([1, 2]), Interval{Float64}([3, 4]))
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
            require(@__MODULE__, :LazySets; fun_name="difference")

            N = promote_type(eltype(X), eltype(Y))
            return EmptySet{N}(1)
        elseif flat_left && !flat_right
            return Interval(h, max(X))
        elseif !flat_left && flat_right
            return Interval(min(X), l)
        else
            require(@__MODULE__, :LazySets; fun_name="difference")

            return UnionSet(Interval(min(X), l), Interval(h, max(X)))
        end
    end
end
