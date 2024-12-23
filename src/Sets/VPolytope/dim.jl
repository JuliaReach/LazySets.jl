"""
    dim(P::VPolytope)

### Output

If `P` is empty, the result is ``-1``.
"""
function dim(P::VPolytope)
    return isempty(P.vertices) ? -1 : @inbounds length(P.vertices[1])
end
