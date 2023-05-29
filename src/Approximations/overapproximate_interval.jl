"""
    overapproximate(S::LazySet, ::Type{<:Interval})

Return the overapproximation of a set with an interval.

### Input

- `S`        -- one-dimensional set
- `Interval` -- target type

### Output

An interval.

### Algorithm

We use the `extrema` function.
"""
function overapproximate(S::LazySet, ::Type{<:Interval})
    @assert dim(S) == 1 "cannot overapproximate a $(dim(S))-dimensional set " *
                        "with an `Interval`"
    l, h = extrema(S, 1)
    return Interval(l, h)
end

"""
    overapproximate(cap::Intersection, ::Type{<:Interval})

Return the overapproximation of a lazy intersection with an interval.

### Input

- `cap`      -- one-dimensional intersection
- `Interval` -- target type

### Output

An interval.

### Algorithm

The algorithm recursively overapproximates the two intersected sets with
intervals and then intersects these. (Note that this fails if the sets are
unbounded.)

For convex sets this algorithm is precise.
"""
function overapproximate(cap::Intersection, ::Type{<:Interval})
    @assert dim(cap) == 1 "cannot overapproximate a $(dim(cap))-dimensional " *
                          "intersection with an `Interval`"
    # TODO this does not work for unbounded sets; better define extrema and
    # then copy the default implementation except if result is empty
    X = overapproximate(cap.X, Interval)
    Y = overapproximate(cap.Y, Interval)
    return intersection(X, Y)
end

"""
    overapproximate(cap::IntersectionArray, ::Type{<:Interval})

Return the overapproximation of an intersection array with an interval.

### Input

- `cap`      -- one-dimensional intersection array
- `Interval` -- target type

### Output

An interval.

### Algorithm

The algorithm recursively overapproximates the intersected sets with intervals
and then intersects these. (Note that this fails if the sets are
unbounded.)

For convex sets this algorithm is precise.
"""
function overapproximate(cap::IntersectionArray, ::Type{<:Interval})
    @assert dim(cap) == 1 "cannot overapproximate a $(dim(cap))-dimensional " *
                          "intersection with an `Interval`"
    a = array(cap)
    # TODO this does not work for unbounded sets; better define extrema and
    # then copy the default implementation except if result is empty
    X = overapproximate(a[1], Interval)
    @inbounds for i in 2:length(a)
        Y = overapproximate(a[i], Interval)
        X = intersection(X, Y)
        if isempty(X)
            return X
        end
    end
    return X
end
