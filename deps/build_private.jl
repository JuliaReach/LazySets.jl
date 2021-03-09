# NOTE: This script is only supposed to be executed by our own build script.

import Pkg

# https://docs.travis-ci.com/user/languages/julia/#dependency-management
@static if VERSION < v"1.4"
    Pkg.add([Pkg.PackageSpec(name="Distributions", version=v"0.22"),
             Pkg.PackageSpec(name="Documenter", version=v"0.25"),
             Pkg.PackageSpec(name="GLPK", version=v"0.13"),
             Pkg.PackageSpec(name="GLPKMathProgInterface", version=v"0.5"),
             Pkg.PackageSpec(name="GR", version=v"0.44"),
             Pkg.PackageSpec(name="IntervalMatrices", version=v"0.6"),
             Pkg.PackageSpec(name="ModelingToolkit", version=v"0.8.0"),
             Pkg.PackageSpec(name="Optim", version=v"0.21"),
             Pkg.PackageSpec(name="TaylorModels", version=v"0.3.6")])
end
