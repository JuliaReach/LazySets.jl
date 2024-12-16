function sample(B::EmptySet, nsamples::Int;
                rng::AbstractRNG=GLOBAL_RNG,
                seed::Union{Int,Nothing}=nothing)
    throw(ArgumentError("cannot sample from an empty set"))
end
