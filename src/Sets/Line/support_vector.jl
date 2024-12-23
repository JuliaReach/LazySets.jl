function σ(d::AbstractVector, L::Line)
    if isapproxzero(dot(d, L.d))
        return L.p
    else
        throw(ArgumentError("the support vector is undefined because the " *
                            "line is unbounded in the given direction"))
    end
end
