# standalone test on computing the mixed volume
include("PHCpack.jl")
verbose = 10
polynomials = ["x^3 + 2*x*y - 1;", "x + y - 1;"]
PHCpack.setDoubleSystem(2, polynomials, verbose)
pols = PHCpack.getDoubleSystem(verbose)
println("the polynomials :")
for pol in pols
    println(pol)
end
mv = PHCpack.mixedVolume(true, verbose)
println("the mixed volume : ", mv)
mv = PHCpack.mixedVolume(false, verbose)
println("the mixed volume : ", mv)
