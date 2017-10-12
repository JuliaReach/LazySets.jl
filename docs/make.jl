using Documenter, LazySets

makedocs(
	modules = [LazySets],
	format = :html,
	sitename = "LazySets.jl",
	pages = [
		"Home" => "index.md",
		"Getting started" => "getting_started.md",
		"Common set representations" => [
            "Balls" => "representations/balls.md",
            "Polygons" => "representations/polygons.md",
            "Other" => "representations/other.md",
        ],
		"Common set operations" => "operations/operations.md",
		"Approximations" => [
            "Cartesian decomposition" => "approximations/cartesian_decomposition.md",
            "Fast polygonal approximation" => "approximations/fast_polygonal_approximation.md",
            "Iterative refinement" => "approximations/iterative_refinement.md"
        ],
		"About" => [
			"Contributing" => "about/CONTRIBUTING.md",
			"License" => "about/license.md",
	    ]
	]
)

deploydocs(
	repo = "github.com/mforets/LazySets.jl",
	target = "build",
	osname = "linux",
	julia  = "0.6",
	deps = nothing,
	make = nothing,
)