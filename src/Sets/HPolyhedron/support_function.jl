"""
    ρ(d::AbstractVector{M}, P::HPoly{N};
      solver=default_lp_solver(M, N)) where {M, N}

Evaluate the support function of a polyhedron in constraint representation in a
given direction.

### Input

- `d`      -- direction
- `P`      -- polyhedron in constraint representation
- `solver` -- (optional, default: `default_lp_solver(M, N)`) the backend used to
              solve the linear program

### Output

The evaluation of the support function for the polyhedron.
If a polytope is unbounded in the given direction, we throw an error.
If a polyhedron is unbounded in the given direction, the result is `Inf`.
"""
function ρ(d::AbstractVector{M}, P::HPoly{N};
           solver=default_lp_solver(M, N)) where {M,N}
    lp, unbounded = σ_helper(d, P, solver)
    if unbounded
        if P isa HPolytope
            error("the support function in direction $(d) is undefined " *
                  "because the polytope is unbounded")
        end
        return N(Inf)
    end
    return dot(d, lp.sol)
end
