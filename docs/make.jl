using Documenter, LazySets

makedocs(
	modules = [LazySets],
	format = :html,
	sitename = "LazySets.jl",
	pages = [
		"Home" => "index.md",
		"Getting Started" => "man/getting_started.md",
		"Common Set Representations" => "man/representations.md", 
        #=
        [
            "Balls" => "representations/hyperrectangles.md",
            "Polygons" => "representations/polygons.md",
            "More types" => "representations/more_types.md"
        ],
        =#
		"Common Set Operations" => "man/operations.md",
		"Approximations" => "man/approximations.md",
        #=
        "Approximations" => [
            "Cartesian decomposition" => "approximations/cartesian_decomposition.md",
            "Box approximations" => "approximations/box_approximations.md",
            "Fast polygonal approximation" => "approximations/fast_polygonal_approximation.md",
            "Iterative refinement" => "approximations/iterative_refinement.md"
        ],
        =#
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