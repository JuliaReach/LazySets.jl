eval(quote
    using .ModelingToolkit: gradient,
                            simplify,
                            Num, # variable like, e.g. x[1]
                            Term # term like, eg. x[1] + x[2] == 1

   """
       _vec(vars::NTuple{L, Union{<:Num, <:Vector{Num}}}) where {L}

   Transform a tuple of operations into one vector of operations.

   ### Input

   - `vars` -- tuple where each element is either a variable-like (`Num`) or a
               vector of variables (`Vector{Num}`)

   ### Output

   A vector of `Operation` obtained by concatenating each tuple component.

   ## Examples

   ```julia
   julia> vars = @variables x[1:2] y
   (Num[x₁, x₂], y)

   julia> LazySets._vec(vars)
   3-element Array{Num,1}:
    x₁
    x₂
    y
   ```
   """
    function _vec(vars::NTuple{L, Union{<:Num, <:Vector{Num}}}) where {L}
        return collect(reduce(vcat, vars))
    end

    # case with a single variable
    _vec(vars::Tuple{Num}) = [vars[1]]

    _get_variables(expr::Term) = ModelingToolkit.get_variables(expr)
    _get_variables(expr::Vector{<:Term}) = unique(reduce(vcat, _get_variables(ex) for ex in expr))
end)

eval(load_modeling_toolkit_hyperplane())
eval(load_modeling_toolkit_halfspace())
eval(load_modeling_toolkit_hpolyhedron())
eval(load_modeling_toolkit_hpolytope())
