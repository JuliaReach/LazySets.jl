@validate function an_element(∅::EmptySet)
    throw(ArgumentError("an empty set does not contain any element"))
end
