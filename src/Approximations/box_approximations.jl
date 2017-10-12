"""
    box_approximation(sf)

Overapproximate a set by a box (hyperrectangle). 

INPUT:

``sf`` -- a set

OUTPUT:

``H`` -- a (tight) hyperrectangle

ALGORITHM:

The center of the hyperrectangle is obtained by averaring the support function
the given set in the canonical directions, and the lengths of the sides can be
recovered from the distance among support functions in the same directions.
"""
function box_approximation(sf::LazySet)::Hyperrectangle
    (n, c, r) = box_approximation_helper(sf)
    return Hyperrectangle(c, r)
end


"""
    box_approximation_symmetric(sf)

Overapproximation of a set by a hyperrectangle which contains the origin.

INPUT:

``sf`` -- a set

OUTPUT:

``H`` -- a symmetric interval around the origin which tightly contains the given set

ALGORITHM:

The center of the box is the origin, and the radius is obtained by computing the
maximum value of the support function evaluated at the canonical directions.
"""
function box_approximation_symmetric(sf::LazySet)::Hyperrectangle
    (n, c, r) = box_approximation_helper(sf)
    return Hyperrectangle(zeros(n), abs.(c) + r)
end
# function alias
symmetric_interval_hull = box_approximation_symmetric


"""
    box_approximation_helper(sf)

Common code of box_approximation and box_approximation_symmetric.

INPUT:

``sf`` -- a set

OUTPUT:

``H`` -- a (tight) hyperrectangle

ALGORITHM:

The center of the hyperrectangle is obtained by averaring the support function
the given set in the canonical directions, and the lengths of the sides can be
recovered from the distance among support functions in the same directions.
"""
@inline function box_approximation_helper(sf::LazySet)
    n = dim(sf)
    c = Vector{Float64}(n)
    r = Vector{Float64}(n)
    dplus = zeros(n)
    dminus = zeros(n)
    @inbounds @simd for i in 1:n
        dplus[i] = 1.0
        dminus[i] = -1.0
        htop = ρ(dplus, sf)
        hbottom = -ρ(dminus, sf)
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

INPUT:

``sf`` -- a set

OUTPUT:

``H`` -- a ball in the infinity norm which tightly contains the given set

ALGORITHM:

The center and radius of the box are obtained by evaluating the support function
of the given set along the canonical directions.
"""
function ballinf_approximation(sf::LazySet)::BallInf
    n = dim(sf)
    c = Vector{Float64}(n)
    r = 0.
    dplus(i::Int64) = [zeros(i-1); 1.; zeros(n-i)]
    dminus(i::Int64) = [zeros(i-1); -1.; zeros(n-i)]
    @inbounds for i in 1:n
        htop = ρ(dplus(i), sf)
        hbottom = -ρ(dminus(i), sf)
        c[i] = (htop+hbottom)/2.
        rcur = (htop-hbottom)/2.
        if (rcur > r)
            r = rcur
        end
    end
    return BallInf(c, r)
end

"""
    radius_approximation(sf)

Approximate radius of a given set.

This is an approximation in the infinity norm. The radius of a BallInf of center
c and radius r can be approximated by ‖c‖ + r√n, where n is the dimension of the 
vectorspace.

INPUT:

``sf`` -- set
"""
function radius_approximation(sf::LazySet)::Float64
    b = ballinf_approximation(sf)
    return norm(b.center) + b.radius * sqrt(dim(sf))
end

"""
    diameter_approximation(sf)

Approximate diameter of a given set.

The diameter is bounded by 2*radius. Relies on radius_approximation.

INPUT:

- ``sf`` -- set
"""
function diameter_approximation(sf::LazySet)::Float64
    return 2.*radius_approximation(sf)
end
