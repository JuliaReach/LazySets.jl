function linear_map(MZ::MatrixZonotope, P::SparsePolynomialZonotope)
	m, n = size(MZ) # dimensions of matrices in MZ
	h = size(P.G, 2)  # order of SPZ
	w = length(MZ.Ai) # num of matrices in MZ
	q = size(P.GI, 2)

	if n != length(center(P))
		throw(DimensionMismatch("Incompatible dimensions:" *
								"size(MZ, 2) = $(n), length(center(P)) = $(length(P.c))"))
	end

	c = center(MZ) * center(P) #center

	# matrix of independent generators
	Gi = Matrix{eltype(MZ)}(undef, m, q * (w + 1))
	Gi[:, 1:q] = center(MZ) * genmat_indep(P)

	# compute matrix of dependendent generators
	G = Matrix{eltype(MZ)}(undef, m, h + w + h * w)
	G[:, 1:h] = center(MZ) * genmat_dep(P)

	# loop to populate G and Gi
	@inbounds for (i, Ai) in generators(MZ)
		G[:, h+i] = Ai * center(P)
		G[:, (h+w+(i-1)*h+1):(h+w+i*h)] = Ai * genmat_dep(P)
		Gi[:, (q*i+1):(q*(i+1))] = Ai * genmat_indep(P)
	end

	# compute exponent
	Imat = Matrix{eltype(expmat(P))}(I, w, w)
	Ê₁, Ê₂, idx = merge_id(P.idx, MZ.idx, P.E, Imat)
	Eₗ = Vector{Matrix{eltype(Ê₁)}}(undef, w)
	@inbounds for l in 1:w
		Eₗ[l] = Ê₂[:, l] * ones(h)' + Ê₁ # fill a matrix of col
	end
	E = hcat(Ê₁, Ê₂, Eₗ...)

	return SparsePolynomialZonotope(c, G, Gi, E, idx)
end
