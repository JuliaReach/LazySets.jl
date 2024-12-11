function high(::EmptySet)
    throw(ArgumentError("the upper bound of an empty set is undefined"))
end

function high(::EmptySet, ::Int)
    throw(ArgumentError("the upper bound of an empty set is undefined"))
end
