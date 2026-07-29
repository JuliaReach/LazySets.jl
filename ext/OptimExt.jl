module OptimExt

using LazySets: HalfSpace, Hyperplane, LazySet, Line2D, isconvex, ρ
using Optim: Brent, optimize
import LazySets: _line_search_optim

function _line_search_optim(ℓ, X::LazySet, H::Union{<:HalfSpace,<:Hyperplane,<:Line2D};
                            kwargs...)
    @assert isconvex(X) "the first set in the intersection must be convex"

    options = Dict(kwargs)

    # Initialization
    a, b = H.a, H.b
    m = -ρ(-ℓ, X)  # `m` is a known lower bound for `f` below (Lemma 3 in paper)
    f(λ) = max(ρ(ℓ - λ[1] * a, X) + λ[1] * b, m)

    if haskey(options, :lower)
        lower = pop!(options, :lower)
    else
        if H isa HalfSpace
            lower = 0.0
        elseif (H isa Hyperplane) || (H isa Line2D)
            lower = -1e6 # "big": TODO relate with f(λ)
        end
    end

    if haskey(options, :upper)
        upper = pop!(options, :upper)
    else
        upper = 1e6 # "big": TODO relate with f(λ)
    end

    if haskey(options, :method)
        method = pop!(options, :method)
    else
        method = Brent()
    end

    # Optimization
    sol = optimize(f, lower, upper; method=method, options...)

    # Recover results
    fmin, λmin = sol.minimum, sol.minimizer
    return (fmin, λmin)
end

end  # module
