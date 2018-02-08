#= conversion between set types =#

# from polygon in H-representation to polygon in H-representation
function convert(::Type{HPOLYGON1},
                 P::HPOLYGON2) where {HPOLYGON1<:AbstractHPolygon,
                                      HPOLYGON2<:AbstractHPolygon}
    return HPOLYGON1(P.constraints)
end

# from polygon in H-representation to polytope in H-representation
function convert(::Type{HPolytope}, P::AbstractHPolygon)
    return HPolytope(P.constraints)
end

# from polytope in H-representation to polygon in H-representation
function convert(::Type{HPOLYGON},
                 P::HPolytope) where {HPOLYGON<:AbstractHPolygon}
    if dim(P) != 2
        error("polytope must be 2D for conversion")
    end
    return HPOLYGON(P.constraints)
end
