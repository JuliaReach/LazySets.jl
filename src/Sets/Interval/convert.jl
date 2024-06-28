function convert(::Type{Interval}, X::LazySet)
    return Interval(convert(IA.Interval, X))
end

function convert(::Type{IA.Interval}, X::LazySet)
    @assert dim(X) == 1 "cannot convert a $(dim(X))-dimensional set to an `Interval`"
    @assert isconvextype(typeof(X)) "cannot convert a non-convex set to an `Interval`"

    l, h = extrema(X, 1)
    return IA.interval(l, h)
end

function convert(::Type{Interval}, x::IA.Interval)
    return Interval(x)
end
