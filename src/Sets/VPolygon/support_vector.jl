# heuristic to define the method used to compute the support vector of a polygon
# in vertex representation; if the number of vertices of the polygon is smaller
# than this value, the brute force method is used; otherwise binary search is used
const BINARY_OR_BRUTE_FORCE = 10

"""
    σ(d::AbstractVector, P::VPolygon)

Return a support vector of a polygon in a given direction.

### Input

- `d` -- direction
- `P` -- polygon in vertex representation

### Output

A support vector in the given direction.
If the direction has norm zero, the first vertex is returned.

### Algorithm

This implementation uses a binary search algorithm when the polygon has more
than $BINARY_OR_BRUTE_FORCE vertices and a brute-force search when it has
$BINARY_OR_BRUTE_FORCE or fewer vertices.
The brute-force search compares the projection of each vector along the given
direction and runs in ``O(n)`` where ``n`` is the number of vertices.
The binary search runs in ``O(log n)`` and we follow
[this implementation](http://geomalgorithms.com/a14-_extreme_pts.html#polyMax_2D())
based on an algorithm described in [1].

[1] Joseph O'Rourke, *Computational Geometry in C (2nd Edition)*.
"""
function σ(d::AbstractVector, P::VPolygon)
    @assert !isempty(P.vertices) "the polygon has no vertices"
    return P.vertices[_σ_helper(d, P)]
end

# return the index of a support vector of P along d from the vertices list of P
function _σ_helper(d::AbstractVector, P::VPolygon)
    vlist = P.vertices
    return _σ_helper(d, vlist)
end

function _σ_helper(d::AbstractVector, vlist::Vector{VN}) where {N,VN<:AbstractVector{N}}
    if length(vlist) > BINARY_OR_BRUTE_FORCE
        return _binary_support_vector(d, vlist)
    else
        return _brute_force_support_vector(d, vlist)
    end
end

function _brute_force_support_vector(d::AbstractVector, P::VPolygon)
    return _brute_force_support_vector(d, P.vertices)
end

function _brute_force_support_vector(d::AbstractVector{M},
                                     vlist::Vector{VT}) where {M,T,VT<:AbstractVector{T}}
    max_idx = 1
    @inbounds max_ρ = dot(d, vlist[1])
    @inbounds for i in 2:length(vlist)
        ρ_i = dot(d, vlist[i])
        if ρ_i > max_ρ
            max_idx = i
            max_ρ = ρ_i
        end
    end
    return max_idx
end

function _binary_support_vector(d::AbstractVector, P::VPolygon)
    return _binary_support_vector(d, P.vertices)
end

# checks if the given vector is pointing toward the given direction
@inline function _similar_direction(u::AbstractVector, v::AbstractVector)
    return dot(u, v) > 0
end

function _binary_support_vector(d::AbstractVector,
                                vlist::Vector{VT}) where {T,VT<:AbstractVector{T}}
    m = length(vlist)
    @assert m > 2 "the number of vertices in the binary-search approach " *
                  "should be at least three, but it is $m"

    @inbounds begin
        # add the first vertex at the end again
        push!(vlist, vlist[1])

        # start chain = [1,n+1] with vlist[n+1] = vlist[1]
        a = 1
        b = m + 1
        A = vlist[2] - vlist[1]
        upA = _similar_direction(d, A)
        if (!upA && !isabove(d, vlist[m], vlist[1]))
            # vlist[1] is the maximum
            pop!(vlist)  # remove the extra point added
            return 1
        end
        while true
            # midpoint of [a,b], and 1<c<n+1
            c = round(Int, (a + b) / 2)
            C = vlist[c + 1] - vlist[c]
            upC = _similar_direction(d, C)
            if (!upC && !isabove(d, vlist[c - 1], vlist[c]))
                # vlist[c] is the maximum
                pop!(vlist)  # remove the extra point added
                return c
            end

            # no max yet, so continue with the binary search
            # pick one of the two subchains [a,c] or [c,b]
            if (upA && upC && !isabove(d, vlist[a], vlist[c])) ||
               (!upA && (upC || (!upC && isabove(d, vlist[a], vlist[c]))))
                a = c
                A = C
                upA = upC
            else
                b = c
            end
            if (b <= a + 1)  # the chain is impossibly small
                pop!(vlist)  # remove the extra point added
                error("something went wrong")
            end
        end
    end
end
