import LazySets.Approximations.overapproximate

# Approximation of a 2D centered unit ball in norm 2
# All vertices v should be like this:
# ‖v‖ >= 1 and ‖v‖ <= 1+ɛ
# Where ɛ is the given error bound
b = Ball2([0., 0.], 1.)
ɛ = .01
p = tovrep(overapproximate(b, ɛ))
for v in p.vl
    @test norm(v) >= 1.
    @test norm(v) <= 1.+ɛ
end
