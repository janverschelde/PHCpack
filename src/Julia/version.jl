"""
    function version(verbose::Int=0)::String

returns the version string of PHCpack,
using ccall() to the function exported by the library libPHCpack.

The input parameter verbose sets the verbose level.
"""
function version(verbose::Int=0)::String
    # The library libPHCpack is made with the Library_Auto_Init flag
    # set to true, so the adainit and adafinal should not be called.
    version = [Cint(0) for i=1:30]
    ptr2version = pointer(version, 1)
    nbr = [Cint(0)]
    ptr2nbr = pointer(nbr, 1)
    flt = [Cdouble(0.0)]
    ptr2flt = pointer(flt, 1)
    p = ccall((:_ada_use_c2phc, LIBPHCPACK), Cint,
              (Cint, Ref{Cint}, Ref{Cint}, Ref{Cdouble}, Cint),
               999, ptr2nbr, ptr2version, ptr2flt, Cint(verbose))
    if verbose > 0
        for c in version
            print(Char(c))
        end
        println()
    end
    result = ""
    for c in version
        result = string(result, Char(c))
    end
    return result
end
