function concretize(qm::QuadraticMap)
    return quadratic_map(qm.Q, concretize(qm.X))
end
