"""
    σ(d::AbstractVector{M}, P::HPoly{N};
      solver=default_lp_solver(M, N) where {M, N}

Return a support vector of a polyhedron in constraint representation in a given
direction.

### Input

- `d`      -- direction
- `P`      -- polyhedron in constraint representation
- `solver` -- (optional, default: `default_lp_solver(M, N)`) the backend used to
              solve the linear program

### Output

The support vector in the given direction.
If a polytope is unbounded in the given direction, we throw an error.
If a polyhedron is unbounded in the given direction, the result contains `±Inf`
entries.
"""
function σ(d::AbstractVector{M}, P::HPoly{N};
           solver=default_lp_solver(M, N)) where {M,N}
    lp, unbounded = σ_helper(d, P, solver)
    if unbounded
        if P isa HPolytope
            error("the support vector in direction $(d) is undefined because " *
                  "the polytope is unbounded")
        end
        return _σ_unbounded_lp(d, P, lp)
    else
        return lp.sol
    end
end

# construct the solution from the solver's ray result
function _σ_unbounded_lp(d, P::HPoly{N}, lp) where {N}
    if isnothing(lp)
        ray = d
    elseif has_lp_infeasibility_ray(lp.model)
        ray = lp.sol  # infeasibility ray is stored as the solution
    else
        error("LP solver did not return an infeasibility ray")
    end

    res = Vector{N}(undef, length(ray))
    e = isempty(P.constraints) ? zeros(N, length(ray)) : an_element(P)
    @inbounds for i in eachindex(ray)
        if isapproxzero(ray[i])
            res[i] = e[i]
        elseif ray[i] > zero(N)
            res[i] = N(Inf)
        else
            res[i] = N(-Inf)
        end
    end
    return res
end

function σ_helper(d::AbstractVector, P::HPoly, solver)
    # represent c = -d as a Vector since GLPK does not accept sparse vectors
    # (see #1011)
    c = to_negative_vector(d)

    (A, b) = tosimplehrep(P)
    if length(b) == 0
        unbounded = true
        lp = nothing
    else
        sense = '<'
        l = -Inf
        u = Inf
        lp = linprog(c, A, sense, b, l, u, solver)
        if is_lp_infeasible(lp.status; strict=true)
            throw(ArgumentError("the support vector is undefined because " *
                                "the polyhedron is empty"))
        elseif is_lp_unbounded(lp.status)
            unbounded = true
        elseif is_lp_optimal(lp.status)
            unbounded = false
        else
            error("got unknown LP status $(lp.status)")
        end
    end
    return (lp, unbounded)
end
