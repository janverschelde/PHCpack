"""
    function setDoublePolynomial(idx::Int, nvr::Int, pol::String,
                                 verbose::Int=0)

sets the polynomial defined by the string pol, in double precision,
with a number of variables no more than nvr, at index idx.

The last input parameter verbose is the verbose level.
"""
function setDoublePolynomial(idx::Int, nvr::Int, pol::String,
                             verbose::Int=0)
    dims = [Cint(length(pol)), Cint(nvr), Cint(idx)]
    ptr2dims = pointer(dims, 1)
    poldata = [Cint(c) for c in pol]
    ptr2poldata = pointer(poldata, 1)
    flt = [Cdouble(0.0)]
    ptr2flt = pointer(flt, 1)
    if verbose > 0
        println("in setDoublePolynomial, setting polynomial ", idx)
    end
    p = ccall((:_ada_use_c2phc, LIBPHCPACK), Cint,
               (Cint, Ref{Cint}, Ref{Cint}, Ref{Cdouble}, Cint),
                76, ptr2dims, ptr2poldata, ptr2flt, Cint(verbose))
end
