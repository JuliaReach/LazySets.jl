@validate function translate(U::Universe, v::AbstractVector)
    return translate!(U, v)  # no need to copy
end

@validate function translate!(U::Universe, v::AbstractVector)
    return U
end
