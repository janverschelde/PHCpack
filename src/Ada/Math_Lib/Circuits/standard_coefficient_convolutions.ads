with Standard_Floating_Vectors;
with Standard_Complex_Vectors;

package Standard_Coefficient_Convolutions is

  function Real_Part ( x : Standard_Complex_Vectors.Link_to_Vector )
                     return Standard_Floating_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns the vector of the real parts of the complex vector x.

  function Imag_Part ( x : Standard_Complex_Vectors.Link_to_Vector )
                     return Standard_Floating_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns the vector of the imaginary parts of the complex vector x.

  function Make_Complex
             ( rpx,ipx : Standard_Floating_Vectors.Link_to_Vector )
             return Standard_Complex_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns the vector of complex numbers,
  --   with real and imaginary parts given in the vectors rpx and ipx.

  procedure Multiply
              ( xr,xi,yr,yi : in Standard_Floating_Vectors.Link_to_Vector;
                zr,zi : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Multiplies the coefficients of the vector x with y,
  --   with real parts in xr, yr, and imaginary parts in xi, yi,
  --   and stores the results in the z, with real and imaginary
  --   parts in zr and zi.

  -- REQUIRED :
  --   All vectors have the same range.

end Standard_Coefficient_Convolutions;
