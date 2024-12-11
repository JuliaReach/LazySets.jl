function norm(::EmptySet, ::Real=Inf)
    throw(ArgumentError("the norm of an empty set is undefined"))
end
