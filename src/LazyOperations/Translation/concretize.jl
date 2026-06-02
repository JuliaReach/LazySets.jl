function concretize(tr::Translation)
    return translate(concretize(tr.X), tr.v)
end
