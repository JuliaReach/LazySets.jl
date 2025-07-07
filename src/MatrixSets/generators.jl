"""
	generators(Z::MatrixZonotope)

Return the generator matrices of the matrix zonotope `Z`.
"""
generators(Z::MatrixZonotope) = Z.Ai
