export RotatedHyperrectangle

struct RotatedHyperrectangle{N, MN<:AbstractMatrix{N},
                             HT<:AbstractHyperrectangle{N}} <: AbstractZonotope{N}
    M::MN
    box::HT

    function RotatedHyperrectangle(M::MN, box::HT) where {N,
            MN<:AbstractMatrix{N}, HT<:AbstractHyperrectangle{N}}
        @assert size(M, 2) == dim(box) "a rotated hyperrectangle of dimension " *
            "$(dim(box)) is incompatible with a matrix of dimension $(size(M))"
        return new{N, MN, HT}(M, box)
    end
end

isoperationtype(::Type{RotatedHyperrectangle}) = false
isconvextype(::Type{RotatedHyperrectangle}) = true

function dim(R::RotatedHyperrectangle)
    return size(R.M, 1)
end

function ρ(d::AbstractVector, R::RotatedHyperrectangle)
    return _ρ_linear_map(d, R.M, R.box)
end

function σ(d::AbstractVector, R::RotatedHyperrectangle)
    return _σ_linear_map(d, R.M, R.box)
end

function center(R::RotatedHyperrectangle)
    return R.M * center(R.box)
end

function generators(R::RotatedHyperrectangle)
    return generators_fallback(R)
end

function genmat(R::RotatedHyperrectangle)
    return R.M * genmat(R.box)
end

function linear_map(M::AbstractMatrix, R::RotatedHyperrectangle)
    @assert size(M, 2) == dim(R) "a linear map of size $(size(M)) is " *
        "incompatible with a set of dimension $(dim(R))"
    if size(M, 1) != size(M, 2)  # TODO generalize
        throw(ArgumentError("non-square linear maps for rotated " *
            "hyperrectangles are not supported yet"))
    end

    return RotatedHyperrectangle(M * R.M, R.box)
end

function vertices_list(R::RotatedHyperrectangle)
    return broadcast(v -> R.M * v, vertices_list(R.box))
end

function constraints_list(R::RotatedHyperrectangle)
    return constraints_list(convert(VPolytope, R))
end
