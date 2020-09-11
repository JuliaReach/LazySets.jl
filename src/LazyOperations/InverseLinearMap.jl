import Base: *, ∈, isempty

export InverseLinearMap

struct InverseLinearMap{N<:Real, S<:LazySet{N},
                 NM, MAT<:AbstractMatrix{NM}} <: AbstractAffineMap{N, S}
    M::MAT
    X::S

    # default constructor with dimension match check
    function InverseLinearMap(M::MAT, X::S) where {N<:Real, S<:LazySet{N}, NM,
                                                   MAT<:AbstractMatrix{NM}}
        @assert dim(X) == size(M, 1) "a linear map of size $(size(M)) cannot " *
            "be applied to a set of dimension $(dim(X))"
        @assert isinvertible(M) "the linear map is not invertible"
        return new{N, S, NM, MAT}(M, X)
    end
end

function ρ(d::AbstractVector{N}, lm::InverseLinearMap{N}; kwargs...) where {N<:Real}
    return ρ(transpose(lm.M) \ d, lm.X; kwargs...)
end

function σ(d::AbstractVector{N}, lm::LinearMap{N}) where {N<:Real}
    return lm.M \ σ(transpose(lm.M) \ d, lm.X)
end
