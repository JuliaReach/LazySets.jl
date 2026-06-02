function concretize(ilm::InverseLinearMap)
    return linear_map(inv(ilm.M), concretize(ilm.X))
end
