module HyperplaneModule

using Reexport, Requires

using ..LazySets: AbstractPolyhedron, AbstractLinearMapAlgorithm,
                  _linear_map_hrep, _non_element_halfspace, _normalize_halfspace
using LinearAlgebra: dot
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Arrays: nonzero_indices
using ReachabilityBase.Commutative: @commutative
using ReachabilityBase.Comparison: _isapprox
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Require: require

@reexport import ..API: an_element, constraints_list, dim, isbounded, isempty,
                        isoperationtype, isuniversal, rand, reflect, distance,
                        ∈, project, ρ, σ, translate
@reexport import ..LazySets: constrained_dimensions, is_hyperplanar, normalize,
                             _is_hyperplane, _linear_map_hrep_helper
@reexport import ..Base: convert
@reexport using ..API

export Hyperplane

include("Hyperplane.jl")

include("an_element.jl")
include("constraints_list.jl")
include("dim.jl")
include("isbounded.jl")
include("isempty.jl")
include("isoperationtype.jl")
include("isuniversal.jl")
include("rand.jl")
include("reflect.jl")
include("distance.jl")
include("in.jl")
include("linear_map.jl")
include("project.jl")
include("support_function.jl")
include("support_vector.jl")
include("translate.jl")

include("constrained_dimensions.jl")
include("is_hyperplanar.jl")
include("normalize.jl")

# ============================================
# Functionality that requires Symbolics
# ============================================
function load_symbolics_hyperplane()
    return quote
        using .Symbolics: Symbolic, Num
        using ..LazySets: _get_variables, _vec

        # returns `(true, sexpr)` if expr represents a hyperplane,
        # where sexpr is the simplified expression sexpr := LHS - RHS == 0
        # otherwise returns `(false, expr)`
        function _is_hyperplane(expr::Symbolic)
            got_hyperplane = Symbolics.operation(expr) == ==
            if got_hyperplane
                # simplify to the form a*x + b == 0
                a, b = Symbolics.arguments(expr)
                sexpr = Symbolics.simplify(a - b)
            end
            return got_hyperplane ? (true, sexpr) : (false, expr)
        end

        """
            Hyperplane(expr::Num, vars=_get_variables(expr); [N]::Type{<:Real}=Float64)

        Return the hyperplane given by a symbolic expression.

        ### Input

        - `expr` -- symbolic expression that describes a hyperplane
        - `vars` -- (optional, default: `_get_variables(expr)`), if a vector of
                    variables is given, use those as the ambient variables with respect
                    to which derivations take place; otherwise, use only the variables
                    that appear in the given expression (but be careful because the
                    order may be incorrect; it is advised to always specify `vars`
                    explicitly)
        - `N`    -- (optional, default: `Float64`) the numeric type of the hyperplane

        ### Output

        A `Hyperplane`.

        ### Examples

        ```jldoctest
        julia> using Symbolics

        julia> vars = @variables x y
        2-element Vector{Num}:
         x
         y

        julia> Hyperplane(x == y)
        Hyperplane{Float64, Vector{Float64}}([1.0, -1.0], -0.0)

        julia> vars = @variables x[1:4]
        1-element Vector{Symbolics.Arr{Num, 1}}:
         x[1:4]

        julia> Hyperplane(x[1] == x[2], x)
        Hyperplane{Float64, Vector{Float64}}([1.0, -1.0, 0.0, 0.0], -0.0)
        ```

        ### Algorithm

        It is assumed that the expression is of the form
        `α*x1 + ⋯ + α*xn + γ == β*x1 + ⋯ + β*xn + δ`.
        This expression is transformed, by rearrangement and substitution, into the
        canonical form `a1 * x1 + ⋯ + an * xn == b`. To identify the coefficients, we
        take derivatives with respect to the ambient variables `vars`. Therefore, the
        order in which the variables appear in `vars` affects the final result. Finally,
        the returned set is the hyperplane with normal vector `[a1, …, an]` and
        displacement `b`.
        """
        function Hyperplane(expr::Num, vars::AbstractVector{Num}=_get_variables(expr);
                            N::Type{<:Real}=Float64)
            valid, sexpr = _is_hyperplane(Symbolics.value(expr))
            if !valid
                throw(ArgumentError("expected an expression of the form `ax == b`, got $expr"))
            end

            # compute the linear coefficients by taking first order derivatives
            coeffs = [N(α.val) for α in Symbolics.gradient(sexpr, collect(vars))]

            # get the constant term by expression substitution
            zeroed_vars = Dict(v => zero(N) for v in vars)
            β = -N(Symbolics.substitute(sexpr, zeroed_vars))

            return Hyperplane(coeffs, β)
        end

        Hyperplane(expr::Num, vars; N::Type{<:Real}=Float64) = Hyperplane(expr, _vec(vars); N=N)
    end
end  # quote / load_symbolics_hyperplane()

# =====================================
# Functionality that requires SymEngine
# =====================================

function load_symengine_hyperplane()
    return quote
        using .SymEngine: Basic
        import .SymEngine: free_symbols
        using ..LazySets: _is_linearcombination

        """
            _is_hyperplane(expr::Expr)

        Determine whether the given expression corresponds to a hyperplane.

        ### Input

        - `expr` -- a symbolic expression

        ### Output

        `true` if `expr` corresponds to a half-space or `false` otherwise.

        ### Examples

        ```jldoctest
        julia> using LazySets: _is_hyperplane

        julia> _is_hyperplane(:(x1 = 0))
        true

        julia> _is_hyperplane(:(x1 <= 0))
        false

        julia> _is_hyperplane(:(2*x1 = 4))
        true

        julia> _is_hyperplane(:(6.1 = 5.3*f - 0.1*g))
        true

        julia> _is_hyperplane(:(2*x1^2 = 4))
        false

        julia> _is_hyperplane(:(x1^2 = 4*x2 - x3))
        false

        julia> _is_hyperplane(:(x1 = 4*x2 - x3))
        true
        ```
        """
        function _is_hyperplane(expr::Expr)::Bool

            # check that there are three arguments
            # these are the comparison symbol, the left hand side and the right hand side
            if (length(expr.args) != 2) || !(expr.head == :(=))
                return false
            end

            # convert to symengine expressions
            lhs = convert(Basic, expr.args[1])

            if :args in fieldnames(typeof(expr.args[2]))
                # treats the 4 in :(2*x1 = 4)
                rhs = convert(Basic, expr.args[2].args[2])
            else
                rhs = convert(Basic, expr.args[2])
            end

            # check if the expression defines a hyperplane
            return _is_linearcombination(lhs) && _is_linearcombination(rhs)
        end

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
            @assert _is_hyperplane(expr) "the expression :(expr) does not correspond to a Hyperplane"

            # get sides of the inequality
            lhs = convert(Basic, expr.args[1])

            # treats the 4 in :(2*x1 = 4)
            rhs = :args in fieldnames(typeof(expr.args[2])) ? convert(Basic, expr.args[2].args[2]) :
                  convert(Basic, expr.args[2])

            # a1 x1 + ... + an xn + K = 0
            eq = lhs - rhs
            if isempty(vars)
                vars = free_symbols(eq)
            end
            K = SymEngine.subs(eq, [vi => zero(N) for vi in vars]...)
            a = convert(Basic, eq - K)

            # convert to numeric types
            K = convert(N, K)
            a = convert(Vector{N}, diff.(a, vars))

            return Hyperplane(a, -K)
        end

        # type-less default Hyperplane conversion
        function convert(::Type{Hyperplane}, expr::Expr; vars::Vector{Basic}=Basic[])
            return convert(Hyperplane{Float64}, expr; vars=vars)
        end

        function free_symbols(expr::Expr, ::Type{<:Hyperplane})
            # get sides of the equality
            lhs = convert(Basic, expr.args[1])

            # treats the 4 in :(2*x1 = 4)
            rhs = :args in fieldnames(typeof(expr.args[2])) ? convert(Basic, expr.args[2].args[2]) :
                  convert(Basic, expr.args[2])
            return free_symbols(lhs - rhs)
        end
    end
end  # quote / load_symengine_hyperplane()

include("init.jl")

end  # module
