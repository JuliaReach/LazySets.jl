"""
    HPolyhedron{N, VN<:AbstractVector{N}} <: AbstractPolyhedron{N}

Type that represents a convex polyhedron in constraint representation, that is,
a finite intersection of half-spaces,
```math
P = ⋂_{i = 1}^m H_i,
```
where each ``H_i = \\{x ∈ ℝ^n : a_i^T x ≤ b_i \\}`` is a
half-space, ``a_i ∈ ℝ^n`` is the normal vector of the ``i``-th
half-space and ``b_i`` is the displacement. The set ``P`` may or may not be
bounded (see also [`HPolytope`](@ref), which assumes boundedness).

### Fields

- `constraints` -- vector of linear constraints
"""
struct HPolyhedron{N,VN<:AbstractVector{N}} <: AbstractPolyhedron{N}
    constraints::Vector{HalfSpace{N,VN}}

    function HPolyhedron(constraints::Vector{HalfSpace{N,VN}}) where {N,VN<:AbstractVector{N}}
        return new{N,VN}(constraints)
    end
end

# constructor with no constraints
function HPolyhedron{N,VN}() where {N,VN<:AbstractVector{N}}
    return HPolyhedron(Vector{HalfSpace{N,VN}}())
end

# constructor with no constraints, given only the numeric type
function HPolyhedron{N}() where {N}
    return HPolyhedron(Vector{HalfSpace{N,Vector{N}}}())
end

# constructor without explicit numeric type, defaults to Float64
function HPolyhedron()
    return HPolyhedron{Float64}()
end

# constructor with constraints of mixed type
function HPolyhedron(constraints::Vector{<:HalfSpace})
    return HPolyhedron(_normal_Vector(constraints))
end

# constructor from a simple constraint representation
function HPolyhedron(A::AbstractMatrix, b::AbstractVector)
    return HPolyhedron(constraints_list(A, b))
end

function load_symbolics_hpolyhedron()
    return quote
        using .Symbolics: Num
        using ..LazySets: _get_variables, _vec, _is_halfspace, _is_hyperplane

        """
            HPolyhedron(expr::Vector{<:Num}, vars=_get_variables(expr);
                        N::Type{<:Real}=Float64)

        Return a polyhedron in constraint representation given by a list of symbolic
        expressions.

        ### Input

        - `expr` -- vector of symbolic expressions that describes each half-space
        - `vars` -- (optional, default: `_get_variables(expr)`), if an array of
                    variables is given, use those as the ambient variables in the set
                    with respect to which derivations take place; otherwise, use only
                    the variables that appear in the given expression (but be careful
                    because the order may be incorrect; it is advised to always pass
                    `vars` explicitly)
        - `N`    -- (optional, default: `Float64`) the numeric type of the returned set

        ### Output

        An `HPolyhedron`.

        ### Examples

        ```jldoctest
        julia> using Symbolics

        julia> vars = @variables x y
        2-element Vector{Num}:
         x
         y

        julia> HPolyhedron([x + y <= 1, x + y >= -1], vars)
        HPolyhedron{Float64, Vector{Float64}}(HalfSpace{Float64, Vector{Float64}}[HalfSpace{Float64, Vector{Float64}}([1.0, 1.0], 1.0), HalfSpace{Float64, Vector{Float64}}([-1.0, -1.0], 1.0)])

        julia> X = HPolyhedron([x == 0, y <= 0], vars)
        HPolyhedron{Float64, Vector{Float64}}(HalfSpace{Float64, Vector{Float64}}[HalfSpace{Float64, Vector{Float64}}([1.0, 0.0], -0.0), HalfSpace{Float64, Vector{Float64}}([-1.0, -0.0], 0.0), HalfSpace{Float64, Vector{Float64}}([0.0, 1.0], -0.0)])
        ```
        """
        function HPolyhedron(expr::Vector{<:Num}, vars::AbstractVector{Num};
                             N::Type{<:Real}=Float64)
            clist = Vector{HalfSpace{N,Vector{N}}}()
            sizehint!(clist, length(expr))
            got_hyperplane = false
            got_halfspace = false
            zeroed_vars = Dict(v => zero(N) for v in vars)
            vars_list = collect(vars)
            for ex in expr
                exval = Symbolics.value(ex)
                got_hyperplane, sexpr = _is_hyperplane(exval)
                if !got_hyperplane
                    got_halfspace, sexpr = _is_halfspace(exval)
                    if !got_halfspace
                        throw(ArgumentError("expected an expression describing either " *
                                            "a half-space of a hyperplane, got $expr"))
                    end
                end

                coeffs = [N(α.val) for α in gradient(sexpr, vars_list)]
                β = -N(Symbolics.substitute(sexpr, zeroed_vars))

                push!(clist, HalfSpace(coeffs, β))
                if got_hyperplane
                    push!(clist, HalfSpace(-coeffs, -β))
                end
            end
            return HPolyhedron(clist)
        end

        function HPolyhedron(expr::Vector{<:Num}; N::Type{<:Real}=Float64)
            return HPolyhedron(expr, _get_variables(expr); N=N)
        end
        function HPolyhedron(expr::Vector{<:Num}, vars; N::Type{<:Real}=Float64)
            return HPolyhedron(expr, _vec(vars); N=N)
        end
    end
end  # quote / load_symbolics_hpolyhedron()
