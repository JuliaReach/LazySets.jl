"""
    difference(X::Interval{N}, Y::Interval) where {N}

Compute the set difference between two intervals.

### Input

- `X` -- first interval
- `Y` -- second interval

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
julia> X = Interval(0, 2); Y = Interval(1, 4);

julia> difference(X, Y)
Interval{Float64}([0, 1])
```
"""
function difference(X::Interval{N}, Y::Interval) where {N}
    Z = intersection(X, Y)
    if isempty(Z)
        return X
    else
        L = Interval(min(X), min(Z))
        R = Interval(max(Z), max(X))

        flat_left = isflat(L)
        flat_right = isflat(R)

        if flat_left && flat_right
            require(@__MODULE__, :LazySets; fun_name="difference")

            return EmptySet{N}(1)
        elseif flat_left && !flat_right
            return R
        elseif !flat_left && flat_right
            return L
        else
            return UnionSet(L, R)
        end
    end
end
