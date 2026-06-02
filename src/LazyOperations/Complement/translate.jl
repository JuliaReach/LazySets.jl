@validate function translate(C::Complement, x::AbstractVector)
    return Complement(translate(C.X, x))
end
