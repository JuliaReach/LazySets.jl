# NOTE: This script is only supposed to be executed by our own build script.
# see: https://docs.travis-ci.com/user/languages/julia/#dependency-management

import Pkg

@static if VERSION > v"1.4"
    Pkg.add([Pkg.PackageSpec(name="Javis", version=v"0.5")])
end
