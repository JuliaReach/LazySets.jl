function convert(::Type{EmptySet}, X::LazySet)
    @assert isempty(X) "cannot convert a nonempty set to an `EmptySet`"

    N = eltype(X)
    n = dim(X)
    return EmptySet{N}(n)
end
