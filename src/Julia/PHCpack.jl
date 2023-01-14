"""
The module PHCpack.jl is a collection of functions to libPHCpack,
the shared object to the library defined by PHCpack.

To check if the shared object libPHCpack is present, type

    PHCpack.version()

and, if present, the version string of PHCpack will appear.
"""
module PHCpack

LIBPHCPACK = "../lib/libPHCpack"

export version, LIBPHCPACK
export setDoubleDimension, getDoubleDimension
export setDoublePolynomial, getDoublePolynomial
export setDoubleSystem, getDoubleSystem
export solveDoubleSystem
export getDoubleSolution
export mixedVolume, stableMixedVolume

include("version.jl")
include("setDoubleDimension.jl")
include("getDoubleDimension.jl")
include("setDoublePolynomial.jl")
include("getDoublePolynomial.jl")
include("setDoubleSystem.jl")
include("getDoubleSystem.jl")
include("solveDoubleSystem.jl")
include("getDoubleSolution.jl")
include("mixedVolume.jl")
include("stableMixedVolume.jl")

end
