using Documenter, LazySets

makedocs(
	modules = [LazySets],
	format = :html,
	sitename = "LazySets.jl",
	pages = [
		"Home" => "index.md",
        "Manual" => Any[
		"Getting Started" => "man/getting_started.md",
		"Common Set Representations" => "lib/representations.md", 
		"Common Set Operations" => "lib/operations.md",
		"Approximations" => "lib/approximations.md",
		"About" => "about.md"]
        "Library" => Any[
		"Getting Started" => "man/getting_started.md",
		"Common Set Representations" => "lib/representations.md", 
		"Common Set Operations" => "lib/operations.md",
		"Approximations" => "lib/approximations.md",
		"About" => "about.md"]        
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