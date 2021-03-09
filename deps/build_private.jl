# NOTE: This script is only supposed to be executed by our own build script.

import Pkg

# https://docs.travis-ci.com/user/languages/julia/#dependency-management
@static if VERSION < v"1.4"
    Pkg.add([Pkg.PackageSpec(name="Documenter", version=v"0.25"),
             Pkg.PackageSpec(name="TaylorModels", version=v"0.3.6")])
end
