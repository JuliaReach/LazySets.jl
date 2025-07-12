@validate function intersection(X::Interval, Y::Interval)
    l = max(min(X), min(Y))
    h = min(max(X), max(Y))
    if l > h
        require(@__MODULE__, :LazySets; fun_name="intersection")

        N = promote_type(eltype(X), eltype(Y))
        return EmptySet{N}(1)
    else
        return Interval(l, h)
    end
end
