"""
# Extended help

    ρ(d::AbstractVector{M}, P::HPoly{N};
      solver=default_lp_solver(M, N)) where {M, N}

### Input

- `solver` -- (optional, default: `default_lp_solver(M, N)`) the backend used to
              solve the linear program

### Output

If `P` is unbounded in the given direction, there are two cases:
- If `P` is an `HPolytope`, we throw an error.
- If `P` is an `HPolyedron`, the result is `Inf`.
"""
@validate function ρ(d::AbstractVector{M}, P::HPoly{N}; solver=default_lp_solver(M, N)) where {M,N}
    lp, unbounded = σ_helper(d, P, solver)
    if unbounded
        if P isa HPolytope
            throw(ArgumentError("the support function in direction $(d) is undefined " *
                                "because the polytope is unbounded"))
        end
        return N(Inf)
    end
    return dot(d, lp.sol)
end
