"""
	center(Z::MatrixZonotope)

Return the center matrix of the matrix zonotope `Z`.
"""
center(Z::MatrixZonotope) = Z.A0
