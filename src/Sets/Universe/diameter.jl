@validate function diameter(::Universe, p::Real=Inf)
    throw(ArgumentError("a universe does not have a diameter"))
end
