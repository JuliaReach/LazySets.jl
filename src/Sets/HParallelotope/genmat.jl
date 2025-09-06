"""
# Extended help

    genmat(P::HParallelotope)

### Algorithm

Let ``P`` be a parallelotope with base vertex ``q`` and list of extremal
vertices with respect to ``q`` given by the set ``\\{v_i\\}`` for
``i = 1, …, n``. Then, the ``i``-th generator of ``P``, represented as the
``i``-th column vector ``G[:, i]``, is given by:

```math
    G[:, i] = \\frac{v_i - q}{2}
```
for ``i = 1, …, n``.
"""
function genmat(P::HParallelotope)
    q = base_vertex(P)
    E = extremal_vertices(P)
    return reduce(hcat, E) / 2 .- q / 2
end
