eval(quote
    using .Symbolics: gradient,
                            simplify,
                            Num,  # variable like, e.g. x[1]
                            Term, # term like, eg. x[1] + x[2] == 1
                            Symbolic,
                            operation,
                            arguments

   """
       _vec(vars::NTuple{L, Union{<:Num, <:Vector{Num}}}) where {L}

   Transform a tuple of operations into one vector of operations.

   ### Input

   - `vars` -- tuple where each element is either variable-like (`Num`) or a
               vector of variables (`Vector{Num}`)

   ### Output

   A vector of `Operation` obtained by concatenating each tuple component.

   ## Examples

   ```julia
   julia> vars = @variables x[1:2] y
   (Num[x₁, x₂], y)

   julia> LazySets._vec(vars)
   3-element Vector{Num}:
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

    # reduce for several variables e.g. when vars = @variables x[1:3] t
    _vec(vars::Vector{Any}) = reduce(vcat, vars)
    _vec(vars::Vector{Num}) = vars
    _vec(vars::Vector{Vector{Num}}) = reduce(vcat, vars)

    _get_variables(expr::Num) = convert(Vector{Num}, Symbolics.get_variables(expr))
    _get_variables(expr::Vector{<:Num}) = unique(reduce(vcat, _get_variables(ex) for ex in expr))
end)

eval(load_symbolics_hyperplane())
eval(load_symbolics_halfspace())
eval(load_symbolics_hpolyhedron())
eval(load_symbolics_hpolytope())
