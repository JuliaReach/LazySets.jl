import Pkg

# https://docs.travis-ci.com/user/languages/julia/#dependency-management
@static if VERSION < v"1.4"
    Pkg.add(Pkg.PackageSpec("Distributions", v"0.22"))
    Pkg.add(Pkg.PackageSpec("Documenter", v"0.25"))
    Pkg.add(Pkg.PackageSpec("GLPK", v"0.13"))
    Pkg.add(Pkg.PackageSpec("GLPKMathProgInterface", v"0.5"))
    Pkg.add(Pkg.PackageSpec("GR", v"0.44"))
    Pkg.add(Pkg.PackageSpec("IntervalArithmetic", v"0.15"))
    Pkg.add(Pkg.PackageSpec("IntervalMatrices", v"0.6"))
    Pkg.add(Pkg.PackageSpec("ModelingToolkit", v"0.8.0"))
    Pkg.add(Pkg.PackageSpec("Optim", v"0.21"))
    Pkg.add(Pkg.PackageSpec("Plots", v"0.28"))
    Pkg.add(Pkg.PackageSpec("RecipesBase", v"0.7"))
    Pkg.add(Pkg.PackageSpec("Requires", v"0.5"))
    Pkg.add(Pkg.PackageSpec("TaylorModels", v"0.1"))
end
