function project(U::Universe{N}, block::AbstractVector{Int}; kwargs...) where {N}
    @assert length(block) <= dim(U) "incompatible dimensions $(length(block)) and $(dim(U))"
    @assert all(1 <= v <= dim(U) for v in block) "invalid dimension in index vector"

    return Universe{N}(length(block))
end
