function Base.convert(::Type{Polygon},
                      U::UnionSet{N,<:Union{AbstractPolytope,Polygon},
                                  <:Union{AbstractPolytope,Polygon}};
                      assume_intersecting::Bool=false) where {N}
    @assert dim(U) == 2 "cannot convert a $(dim(U))-dimensional set to a polygon"

    P = first(U)
    if !(P isa Polygon)
        P = convert(VPolygon, P)
    end
    Q = second(U)
    if !(Q isa Polygon)
        Q = convert(VPolygon, Q)
    end
    return _polygon_union(P, Q; assume_intersecting)
end

# this implementation requires that the sets are given in an order such that the
# next set intersects with the union of the previous sets
function Base.convert(::Type{Polygon},
                      U::UnionSetArray{N,<:Union{AbstractPolytope,Polygon}};
                      assume_intersecting::Bool=false) where {N}
    P = first(U)
    if !(P isa Polygon)
        P = convert(VPolygon, P)
    end

    for Q in U[2:end]
        if !(Q isa Polygon)
            Q = convert(VPolygon, Q)
        end

        P = _polygon_union(P, Q; assume_intersecting)
    end
    return P
end

function _polygon_union(P, Q; assume_intersecting::Bool)
    if !assume_intersecting
        @assert !isdisjoint(P, Q) "the sets must not be disjoint"
    end

    vP = P.vertices
    vQ = Q.vertices

    N = promote_type(eltype(P), eltype(Q))
    vlist = Vector{N}[]

    @inbounds for u in vP
        if u ∉ Q
            push!(vlist, u)
        end
    end
    @inbounds for v in vQ
        if v ∉ P
            push!(vlist, v)
        end
    end
    # TODO very naive iteration over all pairs of line segments
    @inbounds for i in eachindex(vP)
        u1 = vP[i]
        u2 = (i == length(vP)) ? vP[1] : vP[i + 1]
        if u1 == u2
            continue
        end
        for j in eachindex(vQ)
            v1 = vQ[j]
            v2 = (j == length(vQ)) ? vQ[1] : vQ[j + 1]
            if v1 == v2
                continue
            end
            X = intersection(LineSegment(u1, u2), LineSegment(v1, v2))
            if X isa Singleton
                push!(vlist, element(X))
            end
        end
    end

    _sort_ccw!(vlist)

    return Polygon(vlist)
end

function _sort_ccw!(vlist)
    centroid = _centroid_2D(vlist)
    sort!(vlist; by=v -> _angle_from_centroid(v, centroid))
    return vlist
end

function _centroid_2D(vlist)
    @assert !isempty(vlist) "at least one point is required"

    N = eltype(first(vlist))
    x = zero(N)
    y = zero(N)
    @inbounds for v in vlist
        x += v[1]
        y += v[2]
    end
    return [x / length(vlist), y / length(vlist)]
end

function _angle_from_centroid(point, centroid)
    return @inbounds(atan(point[2] - centroid[2], point[1] - centroid[1]))
end
