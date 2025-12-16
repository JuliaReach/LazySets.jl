"""
# Extended help

    in(x::AbstractVector{N}, P::VPolytope{N};
      solver=default_lp_solver(N)) where {N}

### Input

- `solver` -- (optional, default: `default_lp_solver(N)`) the backend used to
              solve the linear program

### Algorithm

We check, using linear programming, the definition of a convex polytope that a
point is in the set if and only if it is a convex combination of the vertices.

Let ``\\{v_j\\}`` be the ``m`` vertices of `P`.
Then we solve the following ``m``-dimensional linear program.

```math
\\max 0 \\text{ s.t. }
⋀_{i=1}^n ∑_{j=1}^m λ_j v_j[i] = x[i]
∧ ∑_{j=1}^m λ_j = 1
∧ ⋀_{j=1}^m λ_j ≥ 0
```
"""
@validate function in(x::AbstractVector{N}, P::VPolytope{N};
                      solver=default_lp_solver(N)) where {N}
    vertices = P.vertices
    m = length(vertices)

    # @validate ensures `m > 0`

    if m == 1
        return _isapprox(x, @inbounds vertices[1])
    end

    n = length(x)
    A = Matrix{N}(undef, n + 1, m)
    for (j, v_j) in enumerate(vertices)
        # ⋀_i Σ_j λ_j v_j[i] = x[i]
        for i in 1:n
            A[i, j] = v_j[i]
        end
        # Σ_j λ_j = 1
        A[n + 1, j] = one(N)
    end
    b = [x; one(N)]
    lbounds = zero(N)
    ubounds = N(Inf)
    sense = '='
    obj = zeros(N, m)
    lp = linprog(obj, A, sense, b, lbounds, ubounds, solver)
    if is_lp_optimal(lp.status)
        return true
    elseif is_lp_infeasible(lp.status)
        return false
    end
    return throw(ArgumentError("unexpected LP solver status: $(lp.status)"))  # COV_EXCL_LINE
end
