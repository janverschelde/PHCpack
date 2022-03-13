# standalone test on setting the dimension
include("JuPHC.jl")
JuPHC.setDoubleDimension(2, 10)
dim = JuPHC.getDoubleDimension(10)
println("the dimension : ", dim)
