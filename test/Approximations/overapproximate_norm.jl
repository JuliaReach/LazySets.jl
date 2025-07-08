using LazySets, Test
using LazySets.Approximations: overapproximate_norm

let
	for _ in 1:10
		Z = rand(Zonotope; dim = 10, num_generators = 10)
		@test overapproximate_norm(Z, 1) ≥ norm(Z, 1)
	end

	Z = Zonotope([1., -1.], zeros(Float64, 2, 2))
	@test norm(Z, 1) == 2
end
