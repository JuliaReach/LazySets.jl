"""
    volume(P::HParallelotope)

Return the volume of a parallelotope in constraint representation.

### Input

- `P` -- parallelotope in constraint representation

### Output

The volume.

### Algorithm

The volume of an ``n``-dimensional parallelotope `P` is ``2^n · |\\det(G)|``,
where ``G`` is the generator matrix of `P`. This can be seen as follows:
The generator matrix transforms the ``n``-dimensional hypercube ``[0, 1]^n``
to a parallelotope of volume ``|\\det(G)|``. Since the representation of a
parallelotope instead transforms the hypercube ``[-1, 1]^n``, this result has to
be doubled for each dimension.
"""
function volume(P::HParallelotope)
    n = dim(P)
    return 2^n * abs(det(genmat(P)))
end
