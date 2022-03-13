"""
    function getDoubleDimension(verbose::Int=0)::Int

returns the number of polynomials in double precision.

The input parameter verbose is the verbose level.
"""
function getDoubleDimension(verbose::Int=0)
    dimstored = [Cint(0)]
    ptr2dimstored = pointer(dimstored, 1)
    nbr = [Cint(0)]
    ptr2nbr = pointer(nbr, 1)
    flt = [Cdouble(0.0)]
    ptr2flt = pointer(flt, 1)
    if verbose > 0
        println("in getDoubleDimension, retrieving dimension ...")
    end
    p = ccall((:_ada_use_c2phc, LIBPHCPACK), Cint,
               (Cint, Ref{Cint}, Ref{Cint}, Ref{Cdouble}, Cint),
                22, ptr2dimstored, ptr2nbr, ptr2flt, Cint(verbose))
    return Int(dimstored[1])
end
