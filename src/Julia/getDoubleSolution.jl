"""
    function getDoubleSolution(idx::Int, verbose::Int=0)::Tuple{Int, String}

returns the string representation of the solution
at the current position, at index idx.

The other input parameter verbose is the verbose level.
"""
function getDoubleSolution(idx::Int, verbose::Int=0)::Tuple{Int, String}

    Cidx = [Cint(idx)]
    ptr2Cidx = pointer(Cidx, 1)

    Csize = [Cint(0)]
    ptr2Csize = pointer(Csize, 1)

    flt = [Cdouble(0.0)]
    ptr2flt = pointer(flt, 1)

    if verbose > 0
        println("in getDoubleSolution, getting the size of ", idx)
    end
    p = ccall((:_ada_use_c2phc, LIBPHCPACK), Cint,
               (Cint, Ref{Cint}, Ref{Cint}, Ref{Cdouble}, Cint),
                525, ptr2Cidx, ptr2Csize, ptr2flt, Cint(verbose))
    size = Int(Csize[1])
    if verbose > 0
        println("in getDoubleSolution, with size = ", size)
    end

    soldata = [Cint(0) for i=1:size]
    ptr2soldata = pointer(soldata, 1)

    if verbose > 0
        println("in getDoubleSolution, getting solution ", idx)
    end
    p = ccall((:_ada_use_c2phc, LIBPHCPACK), Cint,
               (Cint, Ref{Cint}, Ref{Cint}, Ref{Cdouble}, Cint),
                533, ptr2Csize, ptr2soldata, ptr2flt, Cint(verbose))
    if verbose > 0
        for c in soldata
            print(Char(c))
        end
        println()
    end 
    result = ""
    for c in soldata
        result = string(result, Char(c))
    end

    if verbose > 0
        println("in getDoubleSolution, moving the cursor ", idx)
    end
    p = ccall((:_ada_use_c2phc, LIBPHCPACK), Cint,
               (Cint, Ref{Cint}, Ref{Cint}, Ref{Cdouble}, Cint),
                454, ptr2Cidx, ptr2soldata, ptr2flt, Cint(verbose))

    return (Cidx[1], result)
end
