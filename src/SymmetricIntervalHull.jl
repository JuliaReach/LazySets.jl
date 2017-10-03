
struct SymmetricIntervalHull <: LazySet
    sf::LazySet
end

# dimension of a ball in the infinity norm
function dim(SIH::SymmetricIntervalHull)::Int64
    return dim(SIH.sf)
end

function σ(d::Vector{Float64}, SIH::SymmetricIntervalHull)::Vector{Float64}
    #return B.center .+ unit_step.(d) .* B.radius
    return d
end

export SymmetricIntervalHull
