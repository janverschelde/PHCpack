"""
    function setDoubleDimension(dim::Int, verbose::Int=0)

sets the number of polynomials in double precision,
to the value of the first parameter dim.

The other input parameter verbose is the verbose level.
"""
function setDoubleDimension(dim::Int, verbose::Int=0)
    Cdim = [Cint(dim)]
    ptr2Cdim = pointer(Cdim, 1)
    nbr = [Cint(0)]
    ptr2nbr = pointer(nbr, 1)
    flt = [Cdouble(0.0)]
    ptr2flt = pointer(flt, 1)
    if verbose > 0
        println("in setDoubleDimension, initializing dimension ...")
    end
    p = ccall((:_ada_use_c2phc, LIBPHCPACK), Cint,
               (Cint, Ref{Cint}, Ref{Cint}, Ref{Cdouble}, Cint),
                23, ptr2Cdim, ptr2nbr, ptr2flt, Cint(verbose))
end
