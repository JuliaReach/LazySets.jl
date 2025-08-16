function low(::EmptySet)
    throw(ArgumentError("the lower bound of an empty set is undefined"))
end

@validate function low(∅::EmptySet, i::Int)
    throw(ArgumentError("the lower bound of an empty set is undefined"))
end
