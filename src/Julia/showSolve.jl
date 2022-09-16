# standalone test on solving a system
include("PHCpack.jl")
verbose = 10
polynomials = ["x^3 + 2*x*y - 1;", "x + y - 1;"]
PHCpack.setDoubleSystem(2, polynomials, verbose)
pols = PHCpack.getDoubleSystem(verbose)
println("the polynomials :")
for pol in pols
    println(pol)
end
(nbr, roco) = PHCpack.solveDoubleSystem(verbose)
println("number of solutions : ", nbr)
println("root counts :")
println(roco)
idx = 1
for i=1:nbr
    global idx
    (idx, sol) = PHCpack.getDoubleSolution(idx, 0)
    println("SOLUTION ", i, " :")
    println(sol)
end
