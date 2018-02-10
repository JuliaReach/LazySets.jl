#= conversion between set types =#

"""
    convert(::Type{HPOLYGON1}, P::HPOLYGON2) where
        {HPOLYGON1<:AbstractHPolygon, HPOLYGON2<:AbstractHPolygon}

Convert between polygon types in H-representation.

### Input

- `type` -- target type
- `P`    -- source polygon

### Output

The polygon represented as the target type.
"""
function convert(::Type{HPOLYGON1},
                 P::HPOLYGON2) where {HPOLYGON1<:AbstractHPolygon,
                                      HPOLYGON2<:AbstractHPolygon}
    return HPOLYGON1(P.constraints)
end

"""
    convert(::Type{HPolytope}, P::AbstractHPolygon)

Convert from polygon in H-representation to polytope in H-representation.

### Input

- `type` -- target type
- `P`    -- source polygon

### Output

The polygon represented as 2D polytope.
"""
function convert(::Type{HPolytope}, P::AbstractHPolygon)
    return HPolytope(P.constraints)
end

"""
    convert(::Type{HPOLYGON}, P::HPolytope) where {HPOLYGON<:AbstractHPolygon}

Convert from 2D polytope in H-representation to polygon in H-representation.

### Input

- `type` -- target type
- `P`    -- source polytope (must be 2D)

### Output

The 2D polytope represented as polygon.
"""
function convert(::Type{HPOLYGON},
                 P::HPolytope) where {HPOLYGON<:AbstractHPolygon}
    if dim(P) != 2
        error("polytope must be 2D for conversion")
    end
    return HPOLYGON(P.constraints)
end

"""
    convert(::Type{Zonotope}, H::AbstractHyperrectangle{N}) where {N}

Converts a hyperrectangular set to a zonotope.

### Input

- `Zonotope`
- `H` -- hyperrectangular set

### Output

A zonotope.
"""
function convert(::Type{Zonotope}, H::AbstractHyperrectangle{N}) where {N}
    return Zonotope{N}(center(H), diagm(radius_hyperrectangle(H)))
end
