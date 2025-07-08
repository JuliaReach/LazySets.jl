"""
	overapproximate_norm(Z::AbstractZonotope, p::Real=1)

Compute an upper bound on the ``\\ell_p`` norm over a zonotope `Z`.

### Input
- `Z::AbstractZonotope`-- A zonotope or subtype thereof.
- `p::Real`-- The `ℓ_p` norm to overapproximate.

### Output
- `Float64`: An upper bound on ``\\max_{z ∈ Z} \\|z\\|_p``.

"""
function overapproximate_norm(Z::AbstractZonotope, p::Real = 1)
	if p == 1
		return _overapproximate_l1_norm(Z)
	else
		error("an overapproximation of the norm for this value of p=$p is not implemented")
	end
end

"""
	_overapproximate_l1_norm(Z::AbstractZonotope)

Compute an upper bound on the ``\\ell_1`` norm over a zonotope `Z`.

### Notes

The problem ``\\max_{z ∈ Z} \\|z\\|_1`` is NP-hard in general. This function computes a convex relaxation
using a formulation described in [Jordan2021](@citet), which fits a tight upper envelope to the 
absolute value operator. For a zonotope `Z ∈ \\mathbb{R}^d` with `n` generators, this function
has time complexity ``\\mathcal{O}(n ⋅ d)`` compared to ``\\mathcal{O}(2ⁿ ⋅ d)`` for the exact function.
"""
function _overapproximate_l1_norm(Z::AbstractZonotope)
	H = box_approximation(Z)
	ℓ = low(H)
	u = high(H)

	c = center(Z)
	G = generators(Z)

	S⁺ = ℓ .> 0
	S⁻ = u .< 0
	S = .!(S⁺ .| S⁻)

	w = zeros(length(c))
	w[S⁺] .= 1.0
	w[S⁻] .= -1.0
	@. w[S] = (u[S] + ℓ[S]) / (u[S] - ℓ[S])

	const_term = sum(-2 .* u[S] .* ℓ[S] ./ (u[S] .- ℓ[S]))
	linear_max = dot(w, c) + sum(abs.(G' * w))
	return linear_max + const_term
end
