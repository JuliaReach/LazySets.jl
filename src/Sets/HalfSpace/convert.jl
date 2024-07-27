# convert to concrete Vector representation
function convert(::Type{HalfSpace{N,Vector{N}}},
                 hs::HalfSpace{N,<:AbstractVector{N}}) where {N}
    return HalfSpace(Vector(hs.a), hs.b)
end

function load_SymEngine_convert_HalfSpace()
    return quote
        using .SymEngine: Basic
        using ..LazySets: _parse_linear_expression

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
            @assert _is_halfspace(expr) "the expression $expr does not correspond to a half-space"

            # convert to SymEngine expressions
            linexpr, cmp = _parse_halfspace(expr)

            # check sense of the inequality, assuming < or <= by default (checked before)
            got_geq = cmp ∈ (:(>=), :(>))

            # `a1 x1 + ... + an xn + b [cmp] 0` for [cmp] ∈ {<, <=, >, >=}
            a, b = _parse_linear_expression(linexpr, vars, N)

            return got_geq ? HalfSpace(-a, b) : HalfSpace(a, -b)
        end

        # type-less default half-space conversion
        function convert(::Type{HalfSpace}, expr::Expr; vars::Vector{Basic}=Basic[])
            return convert(HalfSpace{Float64}, expr; vars=vars)
        end
    end
end  # load_SymEngine_convert_HalfSpace
