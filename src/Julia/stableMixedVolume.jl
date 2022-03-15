"""
    function stableMixedVolume(demics::Bool=true,
                               verbose::Int=0)::Tuple{Int, Int}

returns the stable mixed volume of the system defined by setDoubleSystem.

The stable mixed volume of a polynomial system equals the number of
isolated solutions (also counting those solutions with zero coordinates)
for a generic choice of the coefficients (not taking the specific
coefficients of the system into account).

The two input parameters are

    demics is the flag to use dynamic enumeration as implemented
    by the software DEMiCs, otherwise MixedVol is called

    verbose is the verbose level parameter.
"""
function stableMixedVolume(demics::Bool=true,
                           verbose::Int=0)::Tuple{Int, Int}
    mv = [Cint(0)]
    ptr2mv = pointer(mv, 1)
    smv = [Cint(0)]
    ptr2smv = pointer(smv, 1)
    flt = [Cdouble(0.0)]
    ptr2flt = pointer(flt, 1)
    if verbose > 0
        if demics
            println("in stableMixedVolume, calling DEMiCs code ...")
        else
            println("in stableMixedVolume, calling MixedVol code ...")
        end
    end
    if demics
        p = ccall((:_ada_use_c2phc, LIBPHCPACK), Cint,
                   (Cint, Ref{Cint}, Ref{Cint}, Ref{Cdouble}, Cint),
                    844, ptr2mv, ptr2smv, ptr2flt, Cint(verbose))
    else
        p = ccall((:_ada_use_c2phc, LIBPHCPACK), Cint,
                   (Cint, Ref{Cint}, Ref{Cint}, Ref{Cdouble}, Cint),
                    79, ptr2mv, ptr2smv, ptr2flt, Cint(verbose))
    end
    return (mv[1], smv[1])
end
