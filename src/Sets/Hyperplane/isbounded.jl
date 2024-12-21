"""
# Extended help

    isbounded(H::Hyperplane)

### Algorithm

The result is `true` iff `H` is one-dimensional.
"""
function isbounded(H::Hyperplane)
    return dim(H) == 1
end
