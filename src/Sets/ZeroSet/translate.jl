function translate(Z::ZeroSet, v::AbstractVector)
    require(@__MODULE__, :LazySets; fun_name="translate")

    @assert length(v) == dim(Z) "cannot translate a $(dim(Z))-dimensional " *
                                "set by a $(length(v))-dimensional vector"

    return Singleton(v)
end
