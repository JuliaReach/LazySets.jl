module HalfSpaceModule

using Reexport, Requires

using ..LazySets: AbstractPolyhedron, LazySet, AbstractLinearMapAlgorithm
import LinearAlgebra
using LinearAlgebra: dot
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Arrays: ismultiple, nonzero_indices, samedir
using ReachabilityBase.Commutative: @commutative
using ReachabilityBase.Comparison: isapproxzero, _isapprox, _leq
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Require: require

@reexport import ..API: an_element, complement, constraints_list, dim,
                        isbounded, isempty, isoperationtype, isuniversal, rand,
                        distance, ∈, permute, project, ρ, σ, translate
@reexport import ..LazySets: constrained_dimensions, normalize, _is_halfspace,
                             _linear_map_hrep_helper
@reexport import ..Base: convert
@reexport using ..API

export HalfSpace, LinearConstraint,
       halfspace_left, halfspace_right,
       iscomplement

include("HalfSpace.jl")

include("an_element.jl")
include("complement.jl")
include("constraints_list.jl")
include("dim.jl")
include("isbounded.jl")
include("isempty.jl")
include("isoperationtype.jl")
include("isuniversal.jl")
include("rand.jl")
include("distance.jl")
include("in.jl")
include("linear_map.jl")
include("permute.jl")
include("project.jl")
include("support_function.jl")
include("support_vector.jl")
include("translate.jl")

include("constrained_dimensions.jl")
include("halfspace_left.jl")
include("halfspace_right.jl")
include("iscomplement.jl")
include("normalize.jl")

include("convert.jl")

"""
    LinearConstraint

Alias for `HalfSpace`
"""
const LinearConstraint = HalfSpace

"""
    is_tighter_same_dir_2D(c1::HalfSpace,
                           c2::HalfSpace;
                           [strict]::Bool=false)

Check if the first of two two-dimensional constraints with equivalent normal
direction is tighter.

### Input

- `c1`     -- first linear constraint
- `c2`     -- second linear constraint
- `strict` -- (optional; default: `false`) check for strictly tighter
              constraints?

### Output

`true` iff the first constraint is tighter.
"""
function is_tighter_same_dir_2D(c1::HalfSpace,
                                c2::HalfSpace;
                                strict::Bool=false)
    @assert dim(c1) == dim(c2) == 2 "the constraints must be two-dimensional"
    @assert samedir(c1.a, c2.a)[1] "the constraints must have the same " *
                                   "normal direction"

    lt = strict ? (<) : (<=)
    if isapproxzero(c1.a[1])
        @assert isapproxzero(c2.a[1])
        return lt(c1.b, c1.a[2] / c2.a[2] * c2.b)
    end
    return lt(c1.b, c1.a[1] / c2.a[1] * c2.b)
end

# TODO: after #2032, #2041 remove use of this function
_normal_Vector(c::HalfSpace) = HalfSpace(convert(Vector, c.a), c.b)
_normal_Vector(C::Vector{<:HalfSpace}) = [_normal_Vector(c) for c in C]
_normal_Vector(P::LazySet) = _normal_Vector(constraints_list(P))

# ============================================
# Functionality that requires Symbolics
# ============================================

function load_symbolics_halfspace()
    return quote
        using .Symbolics: Symbolic, Num
        using ..LazySets: _get_variables, _vec

        # returns `(true, sexpr)` if expr represents a half-space,
        # where sexpr is the simplified expression sexpr := LHS - RHS <= 0
        # otherwise, returns `(false, expr)`
        function _is_halfspace(expr::Symbolic)
            got_halfspace = true

            # find sense and normalize
            op = Symbolics.operation(expr)
            args = Symbolics.arguments(expr)
            if op in (<=, <)
                a, b = args
                sexpr = Symbolics.simplify(a - b)

            elseif op in (>=, >)
                a, b = args
                sexpr = Symbolics.simplify(b - a)

            elseif (op == |) && (Symbolics.operation(args[1]) == <)
                a, b = Symbolics.arguments(args[1])
                sexpr = Symbolics.simplify(a - b)

            elseif (op == |) && (Symbolics.operation(args[2]) == <)
                a, b = Symbolics.arguments(args[2])
                sexpr = Symbolics.simplify(a - b)

            elseif (op == |) && (Symbolics.operation(args[1]) == >)
                a, b = Symbolics.arguments(args[1])
                sexpr = Symbolics.simplify(b - a)

            elseif (op == |) && (Symbolics.operation(args[2]) == >)
                a, b = Symbolics.arguments(args[2])
                sexpr = Symbolics.simplify(b - a)

            else
                got_halfspace = false
            end

            return got_halfspace ? (true, sexpr) : (false, expr)
        end

        """
            HalfSpace(expr::Num, vars=_get_variables(expr); [N]::Type{<:Real}=Float64)

        Return the half-space given by a symbolic expression.

        ### Input

        - `expr` -- symbolic expression that describes a half-space
        - `vars` -- (optional, default: `get_variables(expr)`) if an array of variables
                    is given, use those as the ambient variables in the set with respect
                    to which derivations take place; otherwise, use only the variables
                    that appear in the given expression (but be careful because the
                    order may be incorrect; it is advised to always pass `vars`
                    explicitly; see the examples below for details)
        - `N`    -- (optional, default: `Float64`) the numeric type of the half-space

        ### Output

        A `HalfSpace`.

        ### Examples

        ```jldoctest halfspace_symbolics
        julia> using Symbolics

        julia> vars = @variables x y
        2-element Vector{Num}:
         x
         y

        julia> HalfSpace(x - y <= 2, vars)
        HalfSpace{Float64, Vector{Float64}}([1.0, -1.0], 2.0)

        julia> HalfSpace(x >= y, vars)
        HalfSpace{Float64, Vector{Float64}}([-1.0, 1.0], -0.0)

        julia> vars = @variables x[1:4]
        1-element Vector{Symbolics.Arr{Num, 1}}:
         x[1:4]

        julia> HalfSpace(x[1] >= x[2], x)
        HalfSpace{Float64, Vector{Float64}}([-1.0, 1.0, 0.0, 0.0], -0.0)
        ```

        Be careful with using the default `vars` value because it may introduce a wrong
        order.

        ```@example halfspace_symbolics  # doctest deactivated due to downgrade of Symbolics
        julia> HalfSpace(2x ≥ 5y - 1) # correct
        HalfSpace{Float64, Vector{Float64}}([-2.0, 5.0], 1.0)

        julia> HalfSpace(2x ≥ 5y - 1, vars) # correct
        HalfSpace{Float64, Vector{Float64}}([-2.0, 5.0], 1.0)

        julia> HalfSpace(y - x ≥ 1) # incorrect
        HalfSpace{Float64, Vector{Float64}}([-1.0, 1.0], -1.0)

        julia> HalfSpace(y - x ≥ 1, vars) # correct
        HalfSpace{Float64, Vector{Float64}}([1.0, -1.0], -1.0)
        julia> nothing  # hide
        ```

        ### Algorithm

        It is assumed that the expression is of the form
        `α*x1 + ⋯ + α*xn + γ CMP β*x1 + ⋯ + β*xn + δ`,
        where `CMP` is one among `<`, `<=`, `≤`, `>`, `>=` or `≥`.
        This expression is transformed, by rearrangement and substitution, into the
        canonical form `a1 * x1 + ⋯ + an * xn ≤ b`. The method used to identify the
        coefficients is to take derivatives with respect to the ambient variables `vars`.
        Therefore, the order in which the variables appear in `vars` affects the final
        result. Note in particular that strict inequalities are relaxed as being
        smaller-or-equal. Finally, the returned set is the half-space with normal vector
        `[a1, …, an]` and displacement `b`.
        """
        function HalfSpace(expr::Num, vars::AbstractVector{Num}; N::Type{<:Real}=Float64)
            valid, sexpr = _is_halfspace(Symbolics.value(expr))
            if !valid
                throw(ArgumentError("expected an expression describing a half-space, got $expr"))
            end

            # compute the linear coefficients by taking first-order derivatives
            coeffs = [N(α.val) for α in Symbolics.gradient(sexpr, collect(vars))]

            # get the constant term by expression substitution
            zeroed_vars = Dict(v => zero(N) for v in vars)
            β = -N(Symbolics.substitute(sexpr, zeroed_vars))

            return HalfSpace(coeffs, β)
        end

        HalfSpace(expr::Num; N::Type{<:Real}=Float64) = HalfSpace(expr, _get_variables(expr); N=N)
        HalfSpace(expr::Num, vars; N::Type{<:Real}=Float64) = HalfSpace(expr, _vec(vars); N=N)
    end
end  # quote / load_symbolics_halfspace()

# =====================================
# Functionality that requires SymEngine
# =====================================

function load_symengine_halfspace()
    return quote
        using .SymEngine: Basic
        import .SymEngine: free_symbols
        using ..LazySets: _is_linearcombination

        """
            _is_halfspace(expr::Expr)

        Determine whether the given expression corresponds to a half-space.

        ### Input

        - `expr` -- a symbolic expression

        ### Output

        `true` if `expr` corresponds to a half-space or `false` otherwise.

        ### Examples

        ```jldoctest
        julia> using LazySets: _is_halfspace

        julia> all(_is_halfspace.([:(x1 <= 0), :(x1 < 0), :(x1 > 0), :(x1 >= 0)]))
        true

        julia> _is_halfspace(:(x1 = 0))
        false

        julia> _is_halfspace(:(2*x1 <= 4))
        true

        julia> _is_halfspace(:(6.1 <= 5.3*f - 0.1*g))
        true

        julia> _is_halfspace(:(2*x1^2 <= 4))
        false

        julia> _is_halfspace(:(x1^2 > 4*x2 - x3))
        false

        julia> _is_halfspace(:(x1 > 4*x2 - x3))
        true
        ```
        """
        function _is_halfspace(expr::Expr)::Bool

            # check that there are three arguments
            # these are the comparison symbol, the left hand side and the right hand side
            if (length(expr.args) != 3) || !(expr.head == :call)
                return false
            end

            # check that this is an inequality
            if !(expr.args[1] in [:(<=), :(<), :(>=), :(>)])
                return false
            end

            # convert to symengine expressions
            lhs, rhs = convert(Basic, expr.args[2]), convert(Basic, expr.args[3])

            # check if the expression defines a half-space
            return _is_linearcombination(lhs) && _is_linearcombination(rhs)
        end

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
                vars = free_symbols(eq)
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

        function free_symbols(expr::Expr, ::Type{<:HalfSpace})
            # get sides of the inequality
            lhs, rhs = convert(Basic, expr.args[2]), convert(Basic, expr.args[3])
            return free_symbols(lhs - rhs)
        end
    end
end  # quote / load_symengine_halfspace()

include("init.jl")

end  # module
