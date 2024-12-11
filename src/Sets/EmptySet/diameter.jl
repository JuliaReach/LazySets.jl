function diameter(::EmptySet, ::Real=Inf)
    throw(ArgumentError("the diameter of an empty set is undefined"))
end
