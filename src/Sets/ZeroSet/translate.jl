@validate function translate(Z::ZeroSet, v::AbstractVector)
    require(@__MODULE__, :LazySets; fun_name="translate")

    return Singleton(v)
end
