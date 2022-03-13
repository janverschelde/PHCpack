# standalone test on solving a system
include("JuPHC.jl")
verbose = 10
polynomials = ["x^3 + 2*x*y - 1;", "x + y - 1;"]
JuPHC.setDoubleSystem(2, polynomials, verbose)
pols = JuPHC.getDoubleSystem(verbose)
println("the polynomials :")
for pol in pols
    println(pol)
end
(nbr, roco) = JuPHC.solveDoubleSystem(verbose)
println("number of solutions : ", nbr)
println("root counts :")
println(roco)
idx = 1
for i=1:nbr
    global idx
    (idx, sol) = JuPHC.getDoubleSolution(idx, 0)
    println("SOLUTION ", i, " :")
    println(sol)
end
