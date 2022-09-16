# standalone test on setting a polynomial
include("PHCpack.jl")
verbose = 10
PHCpack.setDoubleDimension(2, verbose)
dim = PHCpack.getDoubleDimension(verbose)
println("the dimension : ", dim)
PHCpack.setDoublePolynomial(1, 2, "x*y - 1;", verbose)
pol = PHCpack.getDoublePolynomial(1, verbose)
println("the polynomial : ", pol)
