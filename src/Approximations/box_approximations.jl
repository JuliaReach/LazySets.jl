"""
    box_approximation(X)

Overapproximate a set by a box (hyperrectangle). 

### Input

`X` -- a lazy set

### Output

`H` -- a (tight) hyperrectangle

### Algorithm

The center of the hyperrectangle is obtained by averaring the support function
of the given set in the canonical directions, and the lengths of the sides can be
recovered from the distance among support functions in the same directions.
"""
function box_approximation(X::LazySet)::Hyperrectangle
    (n, c, r) = box_approximation_helper(X)
    return Hyperrectangle(c, r)
end

"""
    box_approximation_symmetric(X)

Overapproximation of a set by a hyperrectangle which contains the origin.

### Input

`X` -- a lazy set

### Ouptut

`H` -- a symmetric interval around the origin which tightly contains the given set

### Algorithm

The center of the box is the origin, and the radius is obtained by computing the
maximum value of the support function evaluated at the canonical directions.
"""
function box_approximation_symmetric(X::LazySet)::Hyperrectangle
    (n, c, r) = box_approximation_helper(X)
    return Hyperrectangle(zeros(n), abs.(c) + r)
end
# function alias
symmetric_interval_hull = box_approximation_symmetric

"""
    box_approximation_helper(X)

Common code of box_approximation and box_approximation_symmetric.

### Input

`X` -- a lazy set

### Output

`H` -- a (tight) hyperrectangle

### Algorithm

The center of the hyperrectangle is obtained by averaring the support function
the given set in the canonical directions, and the lengths of the sides can be
recovered from the distance among support functions in the same directions.
"""
@inline function box_approximation_helper(X::LazySet)
    n = dim(X)
    c = Vector{Float64}(n)
    r = Vector{Float64}(n)
    dplus = zeros(n)
    dminus = zeros(n)
    @inbounds @simd for i in 1:n
        dplus[i] = 1.0
        dminus[i] = -1.0
        htop = ρ(dplus, X)
        hbottom = -ρ(dminus, X)
        dplus[i] = 0.0
        dminus[i] = 0.0
        c[i] = (htop+hbottom)/2.
        r[i] = (htop-hbottom)/2.
    end
    return n, c, r
end


"""
    ballinf_approximation(S)

Overapproximation of a set by a ball in the infinity norm.

### Input

`X` -- a lazy set

### Output

`H` -- a ball in the infinity norm which tightly contains the given set

### Algorithm

The center and radius of the box are obtained by evaluating the support function
of the given set along the canonical directions.
"""
function ballinf_approximation(X::LazySet)::BallInf
    n = dim(X)
    c = Vector{Float64}(n)
    r = 0.
    dplus = zeros(n)
    dminus = zeros(n)
    @inbounds for i in 1:n
        dplus[i] = 1.0
        dminus[i] = -1.0
        htop = ρ(dplus, X)
        hbottom = -ρ(dminus, X)
        dplus[i] = 0.0
        dminus[i] = 0.0
        c[i] = (htop+hbottom)/2.
        rcur = (htop-hbottom)/2.
        if (rcur > r)
            r = rcur
        end
    end
    return BallInf(c, r)
end

"""
    radius_approximation(X)

Approximate radius of a given set.

### Input

`X` -- a lazy set

### Algorithm

This is an approximation in the infinity norm. The radius of a BallInf of center
c and radius r can be approximated by ‖c‖ + r√n, where n is the dimension of the 
vectorspace.
"""
function radius_approximation(X::LazySet)::Float64
    b = ballinf_approximation(X)
    return norm(b.center) + b.radius * sqrt(dim(X))
end

"""
    diameter_approximation(X)

Approximate diameter of a given set.

### Input

- `X` -- a lazy set

### Algorithm

The diameter is bounded by twice the radius. This function relies on
`radius_approximation`.
"""
function diameter_approximation(X::LazySet)::Float64
    return 2.*radius_approximation(X)
end
