@validate function translate(tr::Translation, x::AbstractVector)
    return Translation(tr.X, tr.v + x)
end
