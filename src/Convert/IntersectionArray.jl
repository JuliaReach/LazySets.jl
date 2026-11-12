"""
    convert(::Type{IntersectionArray}, P::AbstractPolytope)

Convert a two-dimensional polytope to an exact intersection of zonotopes.

### Input

- `IntersectionArray` -- target type
- `P`                 -- two-dimensional polytope

### Output

An `IntersectionArray` whose elements are zonotopes and whose intersection is
equal to `P`.

### Algorithm

We first convert `P` to constraint representation and pair consecutive
constraints. Each pair defines two sides of a containing parallelogram. The
opposite sides are obtained from support-function evaluations in the opposite
directions. If the polygon has an odd number of constraints, the first
constraint is reused in the last pair. Finally, each parallelogram is converted
to a zonotope. See [WanSS18; Theorem 2](@citet).
"""
function convert(::Type{IntersectionArray}, P::AbstractPolytope)
    @assert dim(P) == 2 "the polytope must be two-dimensional"

    if P isa HPolytope && !isbounded(P, false)
        throw(ArgumentError("the polytope must be bounded"))
    end
    Q = convert(HPolygon, P)
    if !isbounded(Q, false)
        throw(ArgumentError("the polytope must be bounded"))
    elseif isempty(Q)
        throw(ArgumentError("the polytope must be nonempty"))
    end

    clist = constraints_list(Q)
    zonotopes = map(1:2:length(clist)) do i
        j = i == length(clist) ? 1 : i + 1
        ci, cj = clist[i], clist[j]

        D = Matrix{eltype(P)}(undef, 2, 2)
        D[1, :] = ci.a
        D[2, :] = cj.a
        if isapproxzero(det(D))
            throw(ArgumentError("cannot pair parallel constraints $i and $j"))
        end

        offset = [ci.b, cj.b, ρ(-ci.a, Q), ρ(-cj.a, Q)]
        H = HParallelotope(D, offset; check_consistency=false)
        return convert(Zonotope, H)
    end
    return IntersectionArray(zonotopes)
end
