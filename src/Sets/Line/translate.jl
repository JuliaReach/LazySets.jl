@validate function translate!(L::Line, v::AbstractVector)
    L.p .+= v
    return L
end
