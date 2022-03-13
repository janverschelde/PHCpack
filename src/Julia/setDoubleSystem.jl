"""
    function setDoubleSystem(nvr::Int, pols::Vector{String},
                             verbose::Int=0)

sets the system defined by the strings in pols, in double precision,
with a number of variables no more than nvr.

The last input parameter verbose is the verbose level.
"""
function setDoubleSystem(nvr::Int, pols::Vector{String},
                         verbose::Int=0)
    dim = length(pols)
    setDoubleDimension(dim, verbose)
    for k=1:dim
        setDoublePolynomial(k, nvr, pols[k], verbose)
    end
end
