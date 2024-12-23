"""
# Extended help

    minkowski_sum(P::VPolygon, Q::VPolygon)

### Algorithm

We treat each edge of the polygons as a vector, attaching them in polar order
(attaching the tail of the next vector to the head of the previous vector). The
resulting polygonal chain will be a polygon, which is the Minkowski sum of the
given polygons. This algorithm assumes that the vertices of `P` and `Q` are
sorted in counter-clockwise fashion and has linear complexity ``O(m+n)``, where
``m`` and ``n`` are the number of vertices of `P` and `Q`, respectively.
"""
function minkowski_sum(P::VPolygon, Q::VPolygon)
    R = _minkowski_sum_vrep_2d(P.vertices, Q.vertices)
    return VPolygon(R)
end
