function load_distributions_sample_from_set()
return quote

using .Distributions: Uniform, Normal

"""
    box_overapproximation(P::LazySet)

Get the values of the box-overapproximation of the set `P`.

### Inputs

- `P`   - lazyset

### Outputs

Matrix with two columns and `n=dims(P)` rows. Each column stands for
one dimension of `P` whereas the first column is the minimum and the second
column the maximum value of the corresponding dimension.
"""
function box_overapproximation(P::LazySet{N}) where {N<:Real}
    dims = dim(P)
    x = Array{N}(undef, dims, 2)
    for j=1:dims
        ej = Arrays.SingleEntryVector(j, dims, 1.0)
        x[j,:] = [-ρ(-ej, P), ρ(ej, P)]
    end
    return x
end


"""
    sample(P::LazySet{N}, nsamples::Int;
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
- `nsamples`    -- number of random samples
- `rng`         -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`        -- (optional, default: `nothing`) seed for reseeding

### Output

A vector of `nsamples` vectors.
"""
function sample(P::LazySet{N}, nsamples::Int=1;
                rng::AbstractRNG=GLOBAL_RNG,
                seed::Union{Int, Nothing}=nothing) where {N<:Real}
    rng = reseed(rng, seed)
    n = dim(P)
    D = Vector{Vector{N}}(undef, nsamples) # preallocate output
    box_ovapp = box_overapproximation(P)
    sampler = [Uniform(box_dim...) for box_dim in eachrow(box_ovapp)]
    @inbounds for i in 1:nsamples
        while true
            w = rand.(Ref(rng), sampler)
            if w ∈ P
                D[i] = w
                break
            end
        end
    end
    return D
end


"""
    Base.rand(P::LazySet;
                   rng::AbstractRNG=GLOBAL_RNG,
                   seed::Union{Int, Nothing}=nothing)

Alias for sample(P::LazySets, 1)
"""
function Base.rand(P::LazySet;
                   rng::AbstractRNG=GLOBAL_RNG,
                   seed::Union{Int, Nothing}=nothing)
    return sample(P, 1)[1]
end

end # quote
end # function load_distributions_sample_from_set()
