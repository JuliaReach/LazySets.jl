function concretize(C::Complement)
    return complement(concretize(C.X))
end
