function concretize(cup::UnionSetArray)
    return UnionSetArray([concretize(X) for X in array(cup)])
end
