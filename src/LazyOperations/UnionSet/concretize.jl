function concretize(cup::UnionSet)
    return UnionSet(concretize(first(cup)), concretize(second(cup)))
end
