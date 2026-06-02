@validate function project(cp::CartesianProduct, block::AbstractVector{Int}; kwargs...)
    n1 = dim(cp.X)
    if block[end] <= n1
        # projection completely in the first block
        return project(cp.X, block; kwargs...)
    elseif block[1] > n1
        # projection completely in the second block
        return project(cp.Y, block .- n1; kwargs...)
    end
    # projection is a new Cartesian product of the block-wise projections
    for (i, bi) in enumerate(block)
        if bi > n1
            X = project(cp.X, block[1:(i - 1)]; kwargs...)
            Y = project(cp.Y, block[i:end] .- n1; kwargs...)
            return CartesianProduct(X, Y)
        end
    end
end

"""
    project(cp::CartesianProduct{N,<:Interval,<:AbstractHyperrectangle},
            block::AbstractVector{Int};
            [kwargs...]) where {N}

Concrete projection of the Cartesian product of an interval and a
hyperrectangular set.

### Input

- `cp`    -- Cartesian product of an interval and a hyperrectangle
- `block` -- block structure, a vector with the dimensions of interest

### Output

A hyperrectangle representing the projection of the Cartesian product `cp` on
the dimensions specified by `block`.
"""
@validate function project(cp::CartesianProduct{N,<:Interval,<:AbstractHyperrectangle},
                           block::AbstractVector{Int}; kwargs...) where {N}
    I = cp.X
    H = cp.Y
    block_vec = collect(block)
    if 1 ∉ block_vec
        block_vec .-= 1
        cH = center(H)
        rH = radius_hyperrectangle(H)
    else
        cH = vcat(center(I), center(H))
        rH = vcat(radius_hyperrectangle(I), radius_hyperrectangle(H))
    end
    return Hyperrectangle(cH[block_vec], rH[block_vec]; check_bounds=false)
end

"""
    project(cp::CartesianProduct{N,<:Interval,<:AbstractZonotope},
            block::AbstractVector{Int};
            [kwargs...]) where {N}

Concrete projection of the Cartesian product of an interval and a zonotopic set.

### Input

- `cp`    -- Cartesian product of an interval and a zonotopic set
- `block` -- block structure, a vector with the dimensions of interest

### Output

A zonotope representing the projection of the Cartesian product `cp` on the
dimensions specified by `block`.
"""
@validate function project(cp::CartesianProduct{N,<:Interval,<:AbstractZonotope},
                           block::AbstractVector{Int}; kwargs...) where {N}
    block_vec = collect(block)
    Z = cp.Y
    if 1 ∉ block_vec
        block_vec .-= 1
    else
        Z = convert(Zonotope, cp)
    end
    M = projection_matrix(block_vec, dim(Z), N)
    return linear_map(M, Z)
end

"""
    project(cp::CartesianProduct{N,<:Interval,<:Union{VPolygon,VPolytope}
            block::AbstractVector{Int};
            [kwargs...]) where {N}

Concrete projection of the Cartesian product of an interval and a set in vertex
representation.

### Input

- `cp`    -- Cartesian product of an interval and a `VPolygon` or a `VPolytope`
- `block` -- block structure, a vector with the dimensions of interest

### Output

A `VPolytope` representing the projection of the Cartesian product `cp` on the
dimensions specified by `block`.
"""
@validate function project(cp::CartesianProduct{N,<:Interval,<:Union{VPolygon,VPolytope}},
                           block::AbstractVector{Int}; kwargs...) where {N}
    I = cp.X
    P = cp.Y
    block_vec = collect(block)
    if 1 ∉ block_vec
        Pout = project(P, block_vec .- 1; kwargs...)
    else
        out = cartesian_product(I, P)
        Pout = project(out, block_vec; kwargs...)
    end
    return Pout
end
