# standalone test on setting a system
include("JuPHC.jl")
verbose = 10
polynomials = ["x^3 + 2*x*y - 1;", "x + y - 1;"]
JuPHC.setDoubleSystem(2, polynomials, verbose)
dim = JuPHC.getDoubleDimension(verbose)
println("the dimension : ", dim)
pols = JuPHC.getDoubleSystem(verbose)
println("the polynomials :")
for pol in pols
    println(pol)
end
