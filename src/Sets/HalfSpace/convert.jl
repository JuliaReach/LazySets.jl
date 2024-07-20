# convert to concrete Vector representation
function convert(::Type{HalfSpace{N,Vector{N}}},
                 hs::HalfSpace{N,<:AbstractVector{N}}) where {N}
    return HalfSpace(Vector(hs.a), hs.b)
end

function load_SymEngine_convert_HalfSpace()
    return quote
        using .SymEngine: Basic

        """
            convert(::Type{HalfSpace{N}}, expr::Expr; vars=nothing) where {N}

        Return a `LazySet.HalfSpace` given a symbolic expression that represents a half-space.

        ### Input

        - `expr` -- a symbolic expression
        - `vars` -- (optional, default: `nothing`): set of variables with respect to which
                    the gradient is taken; if nothing, it takes the free symbols in the given expression

        ### Output

        A `HalfSpace`, in the form `ax <= b`.

        ### Examples

        ```jldoctest convert_halfspace
        julia> convert(HalfSpace, :(x1 <= -0.03))
        HalfSpace{Float64, Vector{Float64}}([1.0], -0.03)

        julia> convert(HalfSpace, :(x1 < -0.03))
        HalfSpace{Float64, Vector{Float64}}([1.0], -0.03)

        julia> convert(HalfSpace, :(x1 > -0.03))
        HalfSpace{Float64, Vector{Float64}}([-1.0], 0.03)

        julia> convert(HalfSpace, :(x1 >= -0.03))
        HalfSpace{Float64, Vector{Float64}}([-1.0], 0.03)

        julia> convert(HalfSpace, :(x1 + x2 <= 2*x4 + 6))
        HalfSpace{Float64, Vector{Float64}}([1.0, 1.0, -2.0], 6.0)
        ```

        You can also specify the set of "ambient" variables, even if not
        all of them appear:

        ```jldoctest convert_halfspace
        julia> using SymEngine: Basic

        julia> convert(HalfSpace, :(x1 + x2 <= 2*x4 + 6), vars=Basic[:x1, :x2, :x3, :x4])
        HalfSpace{Float64, Vector{Float64}}([1.0, 1.0, 0.0, -2.0], 6.0)
        ```
        """
        function convert(::Type{HalfSpace{N}}, expr::Expr; vars::Vector{Basic}=Basic[]) where {N}
            @assert _is_halfspace(expr) "the expression :(expr) does not correspond to a half-space"

            # check sense of the inequality, assuming < or <= by default
            got_geq = expr.args[1] in [:(>=), :(>)]

            # get sides of the inequality
            lhs, rhs = convert(Basic, expr.args[2]), convert(Basic, expr.args[3])

            # a1 x1 + ... + an xn + K [cmp] 0 for cmp in <, <=, >, >=
            eq = lhs - rhs
            if isempty(vars)
                vars = SymEngine.free_symbols(eq)
            end
            K = SymEngine.subs(eq, [vi => zero(N) for vi in vars]...)
            a = convert(Basic, eq - K)

            # convert to numeric types
            K = convert(N, K)
            a = convert(Vector{N}, diff.(a, vars))

            if got_geq
                return HalfSpace(-a, K)
            else
                return HalfSpace(a, -K)
            end
        end

        # type-less default half-space conversion
        function convert(::Type{HalfSpace}, expr::Expr; vars::Vector{Basic}=Basic[])
            return convert(HalfSpace{Float64}, expr; vars=vars)
        end
    end
end  # load_SymEngine_convert_HalfSpace
