"""
    function getDoublePolynomial(idx::Int, verbose::Int=0)

returns the string representation of the polynomial
in double precision, at index idx.

The other input parameter verbose is the verbose level.
"""
function getDoublePolynomial(idx::Int, verbose::Int=0)
    Cidx = [Cint(idx)]
    ptr2Cidx = pointer(Cidx, 1)
    Csize = [Cint(0)]
    ptr2Csize = pointer(Csize, 1)
    flt = [Cdouble(0.0)]
    ptr2flt = pointer(flt, 1)
    if verbose > 0
        println("in getDoublePolynomial, getting the size of ", idx)
    end
    p = ccall((:_ada_use_c2phc, LIBPHCPACK), Cint,
               (Cint, Ref{Cint}, Ref{Cint}, Ref{Cdouble}, Cint),
                600, ptr2Cidx, ptr2Csize, ptr2flt, Cint(verbose))
    size = Int(Csize[1])
    if verbose > 0
        println("in getDoublePolynomial, with size = ", size)
    end
    poldata = [Cint(0) for i=1:size]
    ptr2poldata = pointer(poldata, 1)
    if verbose > 0
        println("in getDoublePolynomial, getting polynomial ", idx)
    end
    p = ccall((:_ada_use_c2phc, LIBPHCPACK), Cint,
               (Cint, Ref{Cint}, Ref{Cint}, Ref{Cdouble}, Cint),
                67, ptr2Cidx, ptr2poldata, ptr2flt, Cint(verbose))
    size = Int(Cidx[1])
    if verbose > 0
        for k = 1:size
            c = poldata[k]
            print(Char(c))
        end
        println()
    end 
    result = ""
    for k = 1:size
        c = poldata[k]
        result = string(result, Char(c))
    end
    return result
end
