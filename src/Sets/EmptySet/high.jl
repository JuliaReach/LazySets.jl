function high(::EmptySet)
    return throw(ArgumentError("the upper bound of an empty set is undefined"))
end

@validate function high(∅::EmptySet, i::Int)
    throw(ArgumentError("the upper bound of an empty set is undefined"))
end
