function low(::EmptySet)
    throw(ArgumentError("the lower bound of an empty set is undefined"))
end

function low(::EmptySet, ::Int)
    throw(ArgumentError("the lower bound of an empty set is undefined"))
end
