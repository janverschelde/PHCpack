# standalone test on setting the dimension
include("PHCpack.jl")
PHCpack.setDoubleDimension(2, 10)
dim = PHCpack.getDoubleDimension(10)
println("the dimension : ", dim)
