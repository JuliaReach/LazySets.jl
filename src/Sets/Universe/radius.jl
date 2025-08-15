@validate function radius(::Universe, p::Real=Inf)
    throw(ArgumentError("a universe does not have a radius"))
end
