"""
This Julia program illustrates the application of ccall()
to call the phctop Ada procedure, which gives access to main.

The library libphcpack is made with the Library_Auto_Init flag
set to true, so the adainit should not be called.
"""

LIBRARY = "../lib/libPHCpack"

# p = ccall((:adainit, LIBRARY), Cvoid, ())
p = ccall((:_ada_phctop, LIBRARY), Cvoid, ())
# p = ccall((:adafinal, LIBRARY), Cvoid, ())
