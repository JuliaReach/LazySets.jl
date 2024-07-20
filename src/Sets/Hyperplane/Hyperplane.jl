"""
    Hyperplane{N, VN<:AbstractVector{N}} <: AbstractPolyhedron{N}

Type that represents a hyperplane of the form ``a⋅x = b``.

### Fields

- `a` -- normal direction (non-zero)
- `b` -- constraint

### Examples

The plane ``y = 0``:

```jldoctest
julia> Hyperplane([0, 1.], 0.)
Hyperplane{Float64, Vector{Float64}}([0.0, 1.0], 0.0)
```
"""
struct Hyperplane{N,VN<:AbstractVector{N}} <: AbstractPolyhedron{N}
    a::VN
    b::N

    function Hyperplane(a::VN, b::N) where {N,VN<:AbstractVector{N}}
        @assert !iszero(a) "a hyperplane needs a non-zero normal vector"
        return new{N,VN}(a, b)
    end
end

function load_Symbolics_Hyperplane()
    return quote
        using .Symbolics: Num
        using ..LazySets: _get_variables, _vec

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
end  # load_Symbolics_Hyperplane
