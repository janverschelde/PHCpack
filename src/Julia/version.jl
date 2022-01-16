"""
This Julia program illustrates the application of ccall()
to call the phctop Ada procedure, which gives access to main.

The library libPHCpack is made with the Library_Auto_Init flag
set to true, so the adainit and adafinal should not be called.
"""

LIBRARY = "../lib/libPHCpack"

version = [Cint(0) for i=1:30]
ptr2version = pointer(version, 1)
nbr = [Cint(0)]
ptr2nbr = pointer(nbr, 1)
flt = [Cdouble(0.0)]
ptr2flt = pointer(flt, 1)
p = ccall((:_ada_use_c2phc, LIBRARY), Cint,
           (Cint, Ref{Cint}, Ref{Cint}, Ref{Cdouble}, Cint),
            999, ptr2nbr, ptr2version, ptr2flt, 1)
# println(version)
for c in version
    print(Char(c))
end
println()
