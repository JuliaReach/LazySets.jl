# NOTE: This script is only supposed to be executed by our own build script.
# see: https://docs.travis-ci.com/user/languages/julia/#dependency-management

import Pkg

@static if VERSION < v"1.4"
    Pkg.add([Pkg.PackageSpec(name="Documenter", version=v"0.26"),
             Pkg.PackageSpec(name="Symbolics", version=v"0.1.2"),
             Pkg.PackageSpec(name="TaylorModels", version=v"0.0.1")])
end

@static if VERSION > v"1.4"
    Pkg.add([Pkg.PackageSpec(name="Symbolics", version=v"0.1.32")])
end
