"""
# Extended help

    volume(B::Ball2)

### Algorithm

This method implements the well-known formula for the volume of an n-dimensional
ball using factorials. For details, see the Wikipedia article
[Volume of an n-ball](https://en.wikipedia.org/wiki/Volume_of_an_n-ball).
"""
function volume(B::Ball2)
    N = eltype(B)
    n = dim(B)
    k = div(n, 2)
    R = B.radius
    if iseven(n)
        vol = N(Base.pi)^k * R^n / factorial(k)
    else
        vol = 2 * factorial(k) * (4 * N(Base.pi))^k * R^n / factorial(n)
    end
    return vol
end
