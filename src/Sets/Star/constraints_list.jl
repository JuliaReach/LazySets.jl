"""
    constraints_list(X::Star)

Return the list of constraints of a star.

### Input

- `X` -- star

### Output

The list of constraints of the star.

### Algorithm

See [`constraints_list(::AbstractAffineMap)`](@ref).
"""
function constraints_list(X::Star)
    require(@__MODULE__, :LazySets; fun_name="constraints_list")

    am = convert(STAR, X)
    return constraints_list(am)
end
