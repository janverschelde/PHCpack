# standalone test on setting a polynomial
include("JuPHC.jl")
verbose = 10
JuPHC.setDoubleDimension(2, verbose)
dim = JuPHC.getDoubleDimension(verbose)
println("the dimension : ", dim)
JuPHC.setDoublePolynomial(1, 2, "x*y - 1;", verbose)
pol = JuPHC.getDoublePolynomial(1, verbose)
println("the polynomial : ", pol)
