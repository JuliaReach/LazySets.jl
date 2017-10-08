using Documenter, LazySets

makedocs(
	modules = [LazySets],
	format = :html,
	sitename = "LazySets.jl",
	pages = [
		"Home" => "index.md",
		"Getting started" => "getting_started.md",
		"Common set representations" => [
			"Unit balls" => "representations/unit_balls.md"
		],
		"Common set operations" => [
			"M-sum" => "operations/msum.md"
		],
		"Approximations" => [
            "decompose" => "approximations/decompose.md"
        ],
		"The lazy_mexp approach" => "approximations/lazy_mexp.md",
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