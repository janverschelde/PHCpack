"""
    function getDoubleSystem(verbose::Int=0)

returns the string representation of the polynomials
in double precision.

The input parameter verbose is the verbose level.
"""
function getDoubleSystem(verbose::Int=0)
    dim = getDoubleDimension(verbose)
    result = []
    for k=1:dim
        pol = getDoublePolynomial(k, verbose)
        push!(result, pol)
    end
    return result
end
