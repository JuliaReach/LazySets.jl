function load_SymEngine_convert_Hyperplane()
    return quote
        using .SymEngine: Basic
        using ..LazySets: _parse_linear_expression

        """
            convert(::Type{Hyperplane{N}}, expr::Expr; vars=nothing) where {N}

        Return a `LazySet.Hyperplane` given a symbolic expression that represents a hyperplane.

        ### Input

        - `expr` -- a symbolic expression
        - `vars` -- (optional, default: `nothing`): set of variables with respect to which
                    the gradient is taken; if nothing, it takes the free symbols in the given expression

        ### Output

        A `Hyperplane`, in the form `ax = b`.

        ### Examples

        ```jldoctest convert_hyperplane
        julia> convert(Hyperplane, :(x1 = -0.03))
        Hyperplane{Float64, Vector{Float64}}([1.0], -0.03)

        julia> convert(Hyperplane, :(x1 + 0.03 = 0))
        Hyperplane{Float64, Vector{Float64}}([1.0], -0.03)

        julia> convert(Hyperplane, :(x1 + x2 = 2*x4 + 6))
        Hyperplane{Float64, Vector{Float64}}([1.0, 1.0, -2.0], 6.0)
        ```

        You can also specify the set of "ambient" variables in the hyperplane, even if not
        all of them appear:

        ```jldoctest convert_hyperplane
        julia> using SymEngine: Basic

        julia> convert(Hyperplane, :(x1 + x2 = 2*x4 + 6), vars=Basic[:x1, :x2, :x3, :x4])
        Hyperplane{Float64, Vector{Float64}}([1.0, 1.0, 0.0, -2.0], 6.0)
        ```
        """
        function convert(::Type{Hyperplane{N}}, expr::Expr; vars::Vector{Basic}=Basic[]) where {N}
            @assert _is_hyperplane(expr) "the expression $expr does not correspond to a Hyperplane"

            # convert to SymEngine expression
            linexpr = _parse_hyperplane(expr)

            # a1 x1 + ... + an xn + b = 0
            a, b = _parse_linear_expression(linexpr, vars, N)

            return Hyperplane(a, -b)
        end

        # type-less default Hyperplane conversion
        function convert(::Type{Hyperplane}, expr::Expr; vars::Vector{Basic}=Basic[])
            return convert(Hyperplane{Float64}, expr; vars=vars)
        end
    end
end  # load_SymEngine_convert_Hyperplane
