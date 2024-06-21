"""
    chebyshev_center_radius(∅::EmptySet; [kwargs]...)

Compute the [Chebyshev center](https://en.wikipedia.org/wiki/Chebyshev_center)
and the corresponding radius of an empty set.

### Input

- `∅`      -- empty set
- `kwargs` -- further keyword arguments (ignored)

### Output

An error.
"""
function chebyshev_center_radius(∅::EmptySet; kwargs...)
    return error("the Chebyshev center and radius of an empty set are undefined")
end
