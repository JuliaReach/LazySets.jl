N = Float64

c = zeros(N, 2)
E1, E2 = diagm(N[-1, 0.5]), N[1 1; 0.5 0.3]
E = [E1, E2];
F2 = N[-0.5 1]';
F = [F2];
G = diagm(N[0.3, 0.3]);
p = PolynomialZonotope(c, E, F, G)

@test dim(p) == 2
@test order(p) == 7//2
@test polynomial_order(p) == 2

# type-specific concrete methods
scale(N(0.5), p)
linear_map(N[1.0 2.0; 2.0 5.0], p)
z = Zonotope(N[1., 2.], N(1.) * eye(N, 2))
minkowski_sum(p, z)
minkowski_sum(z, p)
