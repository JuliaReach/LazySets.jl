using Random, Distributions
using LazySets: dim, GLOBAL_RNG, reseed

function sample(B::Ball2{N}, p::Int=1;
                rng::AbstractRNG=GLOBAL_RNG,
                seed::Union{Int, Nothing}=nothing) where {N}
    
    n = dim(B)
    D = Vector{Vector{N}}(undef, p) # preallocate output
    _sample_unit_nball_muller!(D, n, p, rng=rng)

    # customize for the given ball
    r, c = B.radius, B.center
    @inbounds for i in 1:p
        axpby!(one(N), c, r, D[i])
    end
    return D
end

function _sample_unit_nsphere_muller!(D::Vector{Vector{N}}, n, p;
                                      rng::AbstractRNG=GLOBAL_RNG,
                                      seed::Union{Int, Nothing}=nothing) where {N}
    rng = reseed(rng, seed)
    Zdims = [Normal() for _ in 1:n] # one normal for each dimension
    v = Vector{N}(undef, n) # sample direction
    @inbounds for j in 1:p
        α = zero(N)
        for i in 1:n
            v[i] = rand(rng, Zdims[i])
            α = α + v[i]^2
        end
        D[j] = v ./ sqrt(α)
    end
    return D
end

function _sample_unit_nball_muller!(D::Vector{Vector{N}}, n, p;
                                    rng::AbstractRNG=GLOBAL_RNG,
                                    seed::Union{Int, Nothing}=nothing) where {N}

    rng = reseed(rng, seed)
    Zdims = [Normal() for _ in 1:n] # normals distributions, one for each dimension
    Zrad = Uniform() # distribution to pick random radius
    one_over_n = one(N)/n
    v = Vector{N}(undef, n) # sample direction
    @inbounds for j in 1:p
        α = zero(N)
        for i in 1:n
            v[i] = rand(rng, Zdims[i])
            α = α + v[i]^2
        end
        r = rand(rng, Zrad)
        β = r^one_over_n / sqrt(α)
        D[j] = v .* β
    end
    return D
end
