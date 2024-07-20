"""
    HPolytope{N, VN<:AbstractVector{N}} <: AbstractPolytope{N}

Type that represents a convex polytope in constraint representation, i.e., a
bounded set characterized by a finite intersection of half-spaces,

```math
P = ⋂_{i = 1}^m H_i,
```

where each ``H_i = \\{x ∈ ℝ^n : a_i^T x ≤ b_i \\}`` is a
half-space, ``a_i ∈ ℝ^n`` is the normal vector of the ``i``-th
half-space and ``b_i`` is the displacement.
It is assumed that ``P`` is bounded (see also [`LazySets.HPolyhedron`](@ref),
which does not make such an assumption).

### Fields

- `constraints`       -- vector of linear constraints
- `check_boundedness` -- (optional, default: `false`) flag for checking if the
                         constraints make the polytope bounded; (boundedness is
                         a running assumption for this type)

### Notes

A polytope is a bounded polyhedron.

Boundedness is a running assumption for this type. For performance reasons,
boundedness is not checked in the constructor by default. We also exploit this
assumption, so a boundedness check may not return the answer you would expect.

```jldoctest isbounded
julia> P = HPolytope([HalfSpace([1.0], 1.0)]);  # x <= 1

julia> isbounded(P)  # uses the type assumption and does not actually check
true

julia> isbounded(P, false)  # performs a real boundedness check
false
```
"""
struct HPolytope{N,VN<:AbstractVector{N}} <: AbstractPolytope{N}
    constraints::Vector{HalfSpace{N,VN}}

    function HPolytope(constraints::Vector{HalfSpace{N,VN}};
                       check_boundedness::Bool=false) where {N,VN<:AbstractVector{N}}
        P = new{N,VN}(constraints)
        @assert (!check_boundedness ||
                 isbounded(P, false)) "the polytope is not bounded"
        return P
    end
end

# constructor with no constraints
function HPolytope{N,VN}() where {N,VN<:AbstractVector{N}}
    return HPolytope(Vector{HalfSpace{N,VN}}())
end

# constructor with no constraints, given only the numeric type
function HPolytope{N}() where {N}
    return HPolytope{N,Vector{N}}()
end

# constructor without explicit numeric type, defaults to Float64
function HPolytope()
    return HPolytope{Float64}()
end

# constructor with constraints of mixed type
function HPolytope(constraints::Vector{<:HalfSpace})
    return HPolytope(_normal_Vector(constraints))
end

# constructor from a simple constraint representation
function HPolytope(A::AbstractMatrix, b::AbstractVector;
                   check_boundedness::Bool=false)
    return HPolytope(constraints_list(A, b); check_boundedness=check_boundedness)
end

function load_Polyhedra_HPolytope()
    return quote
        using .Polyhedra: HRep

        """
            HPolytope(P::HRep)

        Return a polytope in constraint representation given an `HRep` polyhedron from
        `Polyhedra.jl`.

        ### Input

        - `P` -- `HRep` polyhedron

        ### Output

        An `HPolytope`.
        """
        function HPolytope(P::HRep)
            return convert(HPolytope, P)
        end
    end
end  # load_Polyhedra_HPolytope

function load_symbolics_hpolytope()
    return quote
        using .Symbolics: Num
        using ..LazySets: _get_variables, _vec

        """
            HPolytope(expr::Vector{<:Num}, vars=_get_variables(expr);
                      [N]::Type{<:Real}=Float64, [check_boundedness]::Bool=false)

        Return the polytope in constraint representation given by a list of symbolic
        expressions.

        ### Input

        - `expr` -- vector of symbolic expressions that describes each constraint
        - `vars` -- (optional, default: `_get_variables(expr)`) if an array of variables
                    is given, use those as the ambient variables in the set with respect
                    to which derivations take place; otherwise, use only the variables
                    that appear in the given expression (but be careful because the
                    order may be incorrect; it is advised to always pass `vars`
                    explicitly)
        - `N`    -- (optional, default: `Float64`) the numeric type of the returned
                    polytope
        - `check_boundedness` -- (optional, default: `false`) flag to check boundedness

        ### Output

        An `HPolytope`.

        ### Examples

        ```jldoctest
        julia> using Symbolics

        julia> vars = @variables x y
        2-element Vector{Num}:
         x
         y

        julia> HPolytope([x <= 1, x >= 0, y <= 1, y >= 0], vars)
        HPolytope{Float64, Vector{Float64}}(HalfSpace{Float64, Vector{Float64}}[HalfSpace{Float64, Vector{Float64}}([1.0, 0.0], 1.0), HalfSpace{Float64, Vector{Float64}}([-1.0, 0.0], 0.0), HalfSpace{Float64, Vector{Float64}}([0.0, 1.0], 1.0), HalfSpace{Float64, Vector{Float64}}([0.0, -1.0], 0.0)])
        ```
        """
        function HPolytope(expr::Vector{<:Num}, vars::AbstractVector{Num};
                           N::Type{<:Real}=Float64, check_boundedness::Bool=false)
            return HPolytope([HalfSpace(ex, vars; N=N) for ex in expr];
                             check_boundedness=check_boundedness)
        end

        function HPolytope(expr::Vector{<:Num}; N::Type{<:Real}=Float64,
                           check_boundedness::Bool=false)
            return HPolytope(expr, _get_variables(expr); N=N,
                             check_boundedness=check_boundedness)
        end

        function HPolytope(expr::Vector{<:Num}, vars; N::Type{<:Real}=Float64,
                           check_boundedness::Bool=false)
            return HPolytope(expr, _vec(vars); N=N, check_boundedness=check_boundedness)
        end
    end
end  # quote / load_symbolics_hpolytope()
