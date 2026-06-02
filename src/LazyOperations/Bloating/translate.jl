@validate function translate(B::Bloating, x::AbstractVector)
    return Bloating(translate(B.X, x), B.ε, B.p)
end
