eval(quote
    using .ModelingToolkit: get_variables,
                            gradient,
                            simplify,
                            Operation

   """
       _vec(vars::NTuple{L, Union{<:Operation, <:Vector{Operation}}}) where {L}

   Transform a tuple of operations into one vector of operations.

   ### Input

   - `vars` -- tuple where each element is either an `Operation` or a vector of operations

   ### Output

   A vector of `Operation` obtained by concatenating each tuple component.

   ## Examples

   ```julia
   julia> vars = @variables x[1:2] y
   (Operation[x₁, x₂], y)

   julia> _vec(vars)
   3-element Array{Operation,1}:
    x₁
    x₂
    y
   ```
   """
    function _vec(vars::NTuple{L, Union{<:Operation, <:Vector{Operation}}}) where {L}
        return collect(reduce(vcat, vars))
    end
end)

eval(load_modeling_toolkit_hyperplane())
eval(load_modeling_toolkit_halfspace())
eval(load_modeling_toolkit_hpolyhedron())
eval(load_modeling_toolkit_hpolytope())
