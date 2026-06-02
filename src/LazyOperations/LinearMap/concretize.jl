function concretize(lm::LinearMap)
    return linear_map(lm.M, concretize(lm.X))
end
