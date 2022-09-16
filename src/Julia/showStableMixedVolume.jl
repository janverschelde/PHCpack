# standalone test on computing the stable mixed volume
include("PHCpack.jl")
verbose = 10
polynomials = ["x^3 + 2*x*y - y^5;", "x + y - x^2;"]
PHCpack.setDoubleSystem(2, polynomials, verbose)
pols = PHCpack.getDoubleSystem(verbose)
println("the polynomials :")
for pol in pols
    println(pol)
end
mv, smv = PHCpack.stableMixedVolume(true, verbose)
println("the mixed volume : ", mv)
println("the stable mixed volume : ", smv)
mv, smv = PHCpack.stableMixedVolume(false, verbose)
println("the mixed volume : ", mv)
println("the stable mixed volume : ", smv)
