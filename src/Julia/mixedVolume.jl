"""
    function mixedVolume(demics::Bool=true, verbose::Int=0)::Int

returns the mixed volume of the system defined by setDoubleSystem.

The mixed volume of a polynomial system equals the number of
isolated solutions (not counting those solutions with zero coordinates)
for a generic choice of the coefficients (not taking the specific
coefficients of the system into account).

The two input parameters are

    demics is the flag to use dynamic enumeration as implemented
    by the software DEMiCs, otherwise MixedVol is called

    verbose is the verbose level parameter.
"""
function mixedVolume(demics::Bool=true, verbose::Int=0)::Int
    mv = [Cint(0)]
    ptr2mv = pointer(mv, 1)
    bpar = [Cint(0)]
    ptr2bpar = pointer(bpar, 1)
    flt = [Cdouble(0.0)]
    ptr2flt = pointer(flt, 1)
    if verbose > 0
        if demics
            println("in mixedVolume, calling DEMiCs code ...")
        else
            println("in mixedVolume, calling MixedVol code ...")
        end
    end
    if demics
        p = ccall((:_ada_use_c2phc, LIBPHCPACK), Cint,
                   (Cint, Ref{Cint}, Ref{Cint}, Ref{Cdouble}, Cint),
                    843, ptr2mv, ptr2bpar, ptr2flt, Cint(verbose))
    else
        p = ccall((:_ada_use_c2phc, LIBPHCPACK), Cint,
                   (Cint, Ref{Cint}, Ref{Cint}, Ref{Cdouble}, Cint),
                    78, ptr2mv, ptr2bpar, ptr2flt, Cint(verbose))
    end
    return mv[1]
end
