@validate function norm(::Universe, p::Real=Inf)
    throw(ArgumentError("a universe does not have a norm"))
end
