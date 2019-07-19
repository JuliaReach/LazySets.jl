 using CDDLib # ..

"""
    minkowski_sum(P::LazySet{N}, Q::LazySet{N};
                  backend=nothing,
                  algorithm=nothing,
                  prune=true) where {N<:Real}

### Input

### Output

### Algorithm

"""
function minkowski_sum(P::LazySet{N}, Q::LazySet{N};
                    backend=nothing,
                    algorithm=nothing,
                    prune=true) where {N<:Real}

    require(:Polyhedra; fun_name="minkowski_sum")
    require(:CDDLib; fun_name="minkowski_sum")

    _minkowski_sum()
end



function minkowski_sum
    if backend == nothing
        if N <: Rational
            backend = CDDLib.Library(:Exact)
         else
             backend = CDDLib.Library()
         end
    end
    if algorithm == nothing
        algorithm = Polyhedra.FourierMotzkin()
    end

    A, b = tosimplehrep(P)
    C, d = tosimplehrep(Q)
    mP, nP = size(A)
    mQ, nQ = size(C)
    E = [zeros(N, mP, nQ) A; C -C]
    f = [b; d]
    PQ = HPolytope(E, f)
    PQ_cdd = polyhedron(PQ, backend=backend)
    W = Polyhedra.eliminate(PQ_cdd, nP+1:2nP, algorithm) |> HPolytope
    if prune
        remove_redundant_constraints!(W) # process bool output?
    end
    return W
end
