@validate function translate(tr::Translation, x::AbstractVector)
    return Translation(translate(tr.X, x))
end
