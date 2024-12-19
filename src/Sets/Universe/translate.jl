function translate(U::Universe, v::AbstractVector)
    return translate!(U, v)  # no need to copy
end

function translate!(U::Universe, v::AbstractVector)
    @assert length(v) == dim(U) "cannot translate a $(dim(U))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    return U
end
