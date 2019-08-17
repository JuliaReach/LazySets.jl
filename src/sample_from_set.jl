function load_distributions_sample_from_set()
return quote

using .Distributions: Uniform, Distribution

"""
    canonical_length(P::LazySet)

Get the canonical length of the box-approximation of the convex set `P`.

### Inputs

- `P`   - lazyset

### Outputs

Matrix with two columns and `n=dims(P)` rows. Each column stands for
one dimension of `P` whereas the first column is the minimum and the second
column the maximum value of the corresponding dimension.
"""
function canonical_length(P::LazySet{N}) where {N<:Real}
    dims = dim(P)
    x = Matrix{N}(undef, dims, 2)
    for j=1:dims
        ej = SingleEntryVector(j, dims, one(N))
        x[j,:] = [-ρ(-ej, P), ρ(ej, P)]
    end
    return x
end

"""
    rand(P::LazySet{N}, n_samples::Int;
           rng::AbstractRNG=GLOBAL_RNG,
           seed::Union{Int, Nothing}=nothing) where {N}

Rejection sampling of an arbitrary LazySet `P` for which the support value
function is defined.

Draw a sample `p` from a uniform distribution of a box-overapproximation of the
original set `P` in all `n` dimensions. The function rejects a drawn sample `p`
and redraws as long as the sample is not contained in the original set `P`,
i.e., `p ∉ P`.

### Input

- `P`           -- lazyset
- `n_samples`    -- number of random samples
- `rng`         -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`        -- (optional, default: `nothing`) seed for reseeding

### Output

A vector of `n_samples` vectors.
"""
function Base.rand(P::LazySet{N}, n_samples::Int;
                rng::AbstractRNG=GLOBAL_RNG,
                seed::Union{Int, Nothing}=nothing) where {N<:Real}
    rng = reseed(rng, seed)
    @assert isbounded(P) "this function requires that the set `P` is bounded, but it is not"

    D = Vector{Vector{N}}(undef, n_samples) # preallocate output
    sample!(D, Sampler(P); rng=rng)
    return D
end

function Base.rand(P::LazySet{N};
                rng::AbstractRNG=GLOBAL_RNG,
                seed::Union{Int, Nothing}=nothing) where {N<:Real}
    return rand(P,1)[1]
end

mutable struct Sampler
    P
    sampler
    distribution
end

function Sampler(P, distribution=Uniform)
    box_lengths = canonical_length(P)
    sampler = [distribution(length...) for length in eachrow(box_lengths)]
    return Sampler(P, sampler, distribution)
end

function sample!(D::Vector{Vector{N}},
               sampler::Sampler;
               rng::AbstractRNG=GLOBAL_RNG) where {N}
    @inbounds for i in 1:length(D)
        while true
            w = rand.(Ref(rng), sampler.sampler)
            if w ∈ sampler.P
                D[i] = w
                break
            end
        end
    end
    nothing
end

end # quote
end # function load_distributions_sample_from_set()
