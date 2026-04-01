function copy(H::Hyperrectangle)
    return Hyperrectangle(copy(H.center), copy(H.radius))
end
