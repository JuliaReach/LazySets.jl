"""
# Extended help

    constraints_list(X::Star)

### Algorithm

See [`constraints_list(::LazySets.AbstractAffineMap)`](@ref).
"""
@validate function constraints_list(X::Star)
    require(@__MODULE__, :LazySets; fun_name="constraints_list")

    am = convert(STAR, X)
    return constraints_list(am)
end
