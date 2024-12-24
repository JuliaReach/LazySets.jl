"""
# Extended help

    dim(P::HPoly)

### Output

If `P` has no constraints, the result is ``-1``.
"""
function dim(P::HPoly)
    return length(P.constraints) == 0 ? -1 : length(P.constraints[1].a)
end
