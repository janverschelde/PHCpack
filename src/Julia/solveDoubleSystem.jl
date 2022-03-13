"""
    function solveDoubleSystem(nbtasks::Int=0, mvfocus::Int=0,
                               verbose::Int=0)::Tuple{Int, String}

solves the system defined by the set operations.

Returns the number of solutions and the root counters string.

The three input parameters are

    nbtasks equals the number of tasks, no multitasking if zero

    mvfocus equals zero by default and all root counts are computed,
    otherwise, the focus is on mixed volumes and polynomial homotopies

    verbose is the verbose parameter.
"""
function solveDoubleSystem(nbtasks::Int=0, mvfocus::Int=0,
                           verbose::Int=0)::Tuple{Int, String}
    pars = [Cint(1), Cint(nbtasks), Cint(mvfocus)]
    ptr2pars = pointer(pars, 1)
    rootcounts = [Cint(0) for k=1:1024]
    ptr2rootcounts = pointer(rootcounts, 1)
    flt = [Cdouble(0.0)]
    ptr2flt = pointer(flt, 1)
    if verbose > 0
        println("in solveDoubleSystem, calling the solver ...")
    end
    p = ccall((:_ada_use_c2phc, LIBPHCPACK), Cint,
               (Cint, Ref{Cint}, Ref{Cint}, Ref{Cdouble}, Cint),
                77, ptr2pars, ptr2rootcounts, ptr2flt, Cint(verbose))
    roco = pars[1]
    nbchars = pars[2]
    result = ""
    for k in 1:nbchars
        result = string(result, Char(rootcounts[k]))
    end
    return (roco, result)
end
