function an_element(::EmptySet)
    return throw(ArgumentError("an empty set does not contain any element"))
end
