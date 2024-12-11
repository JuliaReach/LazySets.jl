function radius(::EmptySet, ::Real=Inf)
    throw(ArgumentError("the radius of an empty set is undefined"))
end
