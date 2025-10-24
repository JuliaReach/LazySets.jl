"""
    constraints_list(P::VPolytope)

### Algorithm

For one- and two-dimensional sets, we respectively convert to an `Interval` or a
`VPolytope` and call the corresponding `constraints_list` function.
For higher-dimensional sets, we use `tohrep` to compute the constraint
representation and call the corresponding `constraints_list` function.
"""
@validate function constraints_list(P::VPolytope)
    require(@__MODULE__, :LazySets; fun_name="constraints_list")

    n = dim(P)
    if n == 1
        return constraints_list(convert(Interval, P))
    elseif n == 2
        return constraints_list(convert(VPolygon, P))
    else
        return constraints_list(tohrep(P))
    end
end
