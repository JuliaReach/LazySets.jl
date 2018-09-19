#= concrete implementations of binary intersections between sets =#

export intersection

"""
    intersection(L1::Line{N}, L2::Line{N}
                )::Union{Singleton{N}, Line{N}, EmptySet{N}} where {N<:Real}

Return the intersection of two 2D lines.

### Input

- `L1` -- first line
- `L2` -- second line

### Output

If the lines are identical, the result is the first line.
If the lines are parallel and not identical, the result is the empty set.
Otherwise the result is the only intersection point.

### Examples

The line ``y = -x + 1`` intersected with the line ``y = x``:

```jldoctest
julia> intersection(Line([-1., 1.], 0.), Line([1., 1.], 1.))
Singleton{Float64}([0.5, 0.5])
julia> intersection(Line([1., 1.], 1.), Line([1., 1.], 1.))
Line{Float64,Array{Float64,1}}([1.0, 1.0], 1.0)

```
"""
function intersection(L1::Line{N}, L2::Line{N}
                     )::Union{Singleton{N}, Line{N}, EmptySet{N}} where {N<:Real}
    b = [L1.b, L2.b]
    a = [transpose(L1.a); transpose(L2.a)]
    try
        # results in LAPACKException or SingularException if parallel
        return Singleton(a \ b)
    catch e
        @assert e isa LAPACKException || e isa SingularException "unexpected " *
            "$(typeof(e)) from LAPACK occurred while intersecting lines:\n$e"
        # lines are parallel
        if an_element(L1) ∈ L2
            # lines are identical
            return L1
        else
            # lines are parallel but not identical
            return EmptySet{N}()
        end
    end
end

"""
    intersection(H1::AbstractHyperrectangle{N},
                 H2::AbstractHyperrectangle{N}
                )::Union{<:Hyperrectangle{N}, EmptySet{N}} where {N<:Real}

Return the intersection of two hyperrectangles.

### Input

- `H1` -- first hyperrectangle
- `H2` -- second hyperrectangle

### Output

If the hyperrectangles do not intersect, the result is the empty set.
Otherwise the result is the hyperrectangle that describes the intersection.

### Algorithm

In each isolated direction `i` we compute the rightmost left border and the
leftmost right border of the hyperrectangles.
If these borders contradict, then the intersection is empty.
Otherwise the result uses these borders in each dimension.
"""
function intersection(H1::AbstractHyperrectangle{N},
                      H2::AbstractHyperrectangle{N}
                     )::Union{Hyperrectangle{N}, EmptySet{N}} where {N<:Real}
    n = dim(H1)
    c1 = center(H1)
    c2 = center(H2)
    r1 = radius_hyperrectangle(H1)
    r2 = radius_hyperrectangle(H2)
    high = Vector{N}(undef, n)
    low = Vector{N}(undef, n)
    for i in 1:n
        high1 = c1[i] + r1[i]
        low1 = c1[i] - r1[i]
        high2 = c2[i] + r2[i]
        low2 = c2[i] - r2[i]
        high[i] = min(high1, high2)
        low[i] = max(low1, low2)
        if high[i] < low[i]
            return EmptySet{N}()
        end
    end
    return Hyperrectangle(high=high, low=low)
end

function load_polyhedra_concrete_intersection() # function to be loaded by Requires
return quote

export intersection

function default_intersection_output_type(P1::Union{HPolytope{N}, VPolytope{N}},
                                          P2::Union{HPolytope{N}, VPolytope{N}}) where {N}
    T = (P1 isa VPolytope && P2 isa VPolytope) ? VPolytope : HPolytope
    return T
end

"""
    intersection(P1::Union{HPolytope{N}, VPolytope{N}},
                 P2::Union{HPolytope{N}, VPolytope{N}},
                 [backend]=default_polyhedra_backend(N),
                 [prunefunc]=removehredundancy!,
                 [output_type]=default_intersection_output_type(P1, P2)) where {N}

Compute the intersection of two polytopes in either H-representation or
V-representation.

### Input

- `P1`          -- polytope
- `P2`          -- another polytope
- `backend`     -- (optional, default: `default_polyhedra_backend(N)`) the polyhedral
                   computations backend, see [Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/latest/installation.html#Getting-Libraries-1)
                   for further information
- `prunefunc`   -- (optional, default: `removehredundancy!`) function to post-process
                    the output of `intersect`
- `output_type` -- (optional, default: `default_intersection_output_type(P1, P2)`)
                   choose the type for the output polytope 

### Output

The polytope obtained by the intersection of `P1` and `P2`.

The type of the output polytope can be passed optionally. By default, if both
`P1` and `P2` are polytopes in V-representation, the output is a polytope in
V-representation. Otherwise, a polytope in H-representation is returned. 
"""
function intersection(P1::Union{HPolytope{N}, VPolytope{N}},
                      P2::Union{HPolytope{N}, VPolytope{N}},
                      backend=default_polyhedra_backend(N),
                      prunefunc=removehredundancy!,
                      output_type=default_intersection_output_type(P1, P2)) where {N}

    P1 = polyhedron(P1, backend)
    P2 = polyhedron(P2, backend)
    Pint = Polyhedra.intersect(P1, P2)
    prunefunc(Pint)
    return output_type(Pint)
end

end # quote
end # function load_polyhedra_concrete_intersection()
