# standalone test on setting a system
include("PHCpack.jl")
verbose = 10
polynomials = ["x^3 + 2*x*y - 1;", "x + y - 1;"]
PHCpack.setDoubleSystem(2, polynomials, verbose)
dim = PHCpack.getDoubleDimension(verbose)
println("the dimension : ", dim)
pols = PHCpack.getDoubleSystem(verbose)
println("the polynomials :")
for pol in pols
    println(pol)
end
