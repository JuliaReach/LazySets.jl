"""
	ngens(Z::MatrixZonotope)

Return the number of generators of the matrix zonotope `Z`.
"""
function ngens(Z::MatrixZonotope)
	return length(generators(Z))
end
