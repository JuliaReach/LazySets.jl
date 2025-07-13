"""
	overapproximate(lm::LinearMap{N, SparsePolynomialZonotope{N}, NM, MAT}) where {N, NM,
	                 MAT <: MatrixZonotope{NM}}

Overapproximate the linear map of a sparse polynomial zonotope through a matrix zonotope,
following Proposition 1 of []@citet.

### Input

- `lm` -- a linear map of a sparse polynomial zonotope through a matrix zonotope

### Output

A sparse polynomial zonotope overapproximating the linear map.

"""
function overapproximate(lm::LinearMap{N,S,NM,MAT}) where {N,S<:SparsePolynomialZonotope{N},NM,
                                                           MAT<:MatrixZonotope{NM}}
	MZ = matrix(lm)
	P = set(lm)
	T = promote_type(eltype(MZ), eltype(P))

	m, n = size(MZ)
	w = ngens(MZ)
	h = ngens_dep(P)
	q = ngens_indep(P)

	if n != dim(P)
		throw(DimensionMismatch("incompatible dimensions:" *
								"size(MZ) = $(size(MZ)), dim(P) = $q"))
	end

	c = center(MZ) * center(P)

	# matrix of independent generators
	Gi = Matrix{T}(undef, m, q * (w + 1))
	Gi[:, 1:q] = center(MZ) * genmat_indep(P)

	# compute matrix of dependendent generators
	G = Matrix{T}(undef, m, h + w + h * w)
	G[:, 1:h] = center(MZ) * genmat_dep(P)

	# loop to populate G and Gi
	@inbounds for (i, A) in enumerate(generators(MZ))
		G[:, h+i] = A * center(P)
		G[:, (h+w+(i-1)*h+1):(h+w+i*h)] = A * genmat_dep(P)
		Gi[:, (q*i+1):(q*(i+1))] = A * genmat_indep(P)
	end

	# compute exponent
	Imat = Matrix{Int}(I, w, w)
	Ê₁, Ê₂, idx = merge_id(indexvector(P), indexvector(MZ), expmat(P), Imat)
	pₖ = size(Ê₁, 1)
	E = Matrix{eltype(Ê₁)}(undef, pₖ, h + w + h * w)
	E[:, 1:h] = Ê₁
	E[:, (h+1):(h+w)] = Ê₂
	@inbounds for l in 1:w
		cstart = (h + w) + (l - 1) * h + 1
		cend = (h + w) + l * h
		E[:, cstart:cend] = Ê₂[:, l] * ones(1, h) .+ Ê₁
	end

	return SparsePolynomialZonotope(c, G, Gi, E, idx)
end


function _taylor_expmap(A::T, B::MatrixZonotope, P::SparsePolynomialZonotope,
                        k::Int) where {T<:Union{MatrixZonotope,Nothing}}
	N = promote_type(eltype(A), eltype(B), eltype(P))

	#  B^i * X / i!
	Bpowers = Vector{SparsePolynomialZonotope{N}}(undef, k)
	f = N(1)
	@inbounds for i in 1:k
		f *= i
		P = overapproximate((1 / f) * B * P)
        P = remove_redundant_generators(P)
		Bpowers[i] = P
	end

	if isnothing(A)
		return reduce(exact_sum, Bpowers)
	end

	# A^i * (B^i * X)/ i!
	result = overapproximate(A * Bpowers[1])
	@inbounds for i in 2:k
		tmp = Bpowers[i]
		for _ in 1:i
			tmp = overapproximate(A * tmp) # compute A^i B^i P
			tmp = remove_redundant_generators(tmp)
		end
		exact_sum(result, tmp)
	end

	return result
end

# helper to pull out A,B
get_factors(MZP::MatrixZonotope) = (nothing, MZP)
get_factors(MZP::MatrixZonotopeProduct) = (MZP.A, MZP.B)

"""
    overapproximate(em::ExponentialMap{N, SparsePolynomialZonotope{N}, NM, MAT},
	                 k::Int = 2) where {N, NM, MAT <: AbstractMatrixZonotope{NM}}

Overapproximate the exponential map of a sparse polynomial zonotope through a product of matrix 
zonotopes, following Proposition 1 of []@citet.

### Input

- `em` -- an expontial map of a sparse polynomial zonotope through a product of matrix zonotopes

### Output

A sparse polynomial zonotope overapproximating the exponential map.

"""
function overapproximate(em::ExponentialMap{N,S,NM,MAT},
                         k::Int=2) where {N,S<:SparsePolynomialZonotope{N},NM,
                                          MAT<:AbstractMatrixZonotope{NM}}
	MZP  = matrix(em)
	A, B = get_factors(MZP)
	P    = set(em)
	n    = size(B, 1)    # A could be nothing, so size(B)

    ABn = (A === nothing ? overapproxi_norm(B, Inf) :
           overapproximate_norm(A, Inf) * overapproximate_norm(B, Inf))
	ϵ = ABn / (k + 2)
	if ϵ > 1
		@warn "κ should be chosen such that ϵ<1 " ϵ
	end

	σ = _taylor_expmap(A, B, P, k)
	ε = IntervalMatrix(fill(IA.interval(-1, 1), n, n))
	ε *= ABn^(k + 1) / (factorial(k + 1) * (1 - ϵ))

	Zp = overapproximate(P)
	rhs = overapproximate(ε * Zp, Zonotope)

	return minkowski_sum(σ, rhs)
end
