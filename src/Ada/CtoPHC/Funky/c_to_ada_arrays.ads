with C_Integer_Arrays;                    use C_Integer_Arrays;
with C_Double_Arrays;                     use C_Double_Arrays;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_Matrices;

package C_to_Ada_Arrays is

-- DESCRIPTION :
--   This package offers conversion operators between arrays of
--   integers and doubles, to interface with C, and the vector types
--   defined in PHCpack.  Note that C arrays always start at 0.

  function Convert ( v : C_Integer_Array ) 
                   return Standard_Integer_Vectors.Vector;
  function Convert ( v : Standard_Integer_Vectors.Vector )
                   return C_Integer_Array;

  -- DESCRIPTION :
  --   Conversion between arrays of C integers and Ada integers.

  function Convert ( v : C_Double_Array )
                   return Standard_Floating_Vectors.Vector;
  function Convert ( v : Standard_Floating_Vectors.Vector )
                   return C_Double_Array;

  -- DESCRIPTION :
  --   Conversion between arrays of C doubles and Ada doubles.

  function Convert ( v : C_Double_Array )
                   return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Converts a double array into an array of complex numbers,
  --   the array on return starts at index 1.

  -- REQUIRED : v'length mod 2 = 0.

  function Convert ( v : Standard_Complex_Vectors.Vector )
                   return C_Double_Array;
  function Convert ( v : DoblDobl_Complex_Vectors.Vector )
                   return C_Double_Array;
  function Convert ( v : QuadDobl_Complex_Vectors.Vector )
                   return C_Double_Array;

  -- DESCRIPTION :
  --   Creates an array of doubles with real and imaginary parts of
  --   the complex numbers in v.

  function Convert ( m : Standard_Complex_Matrices.Matrix )
                   return C_Double_Array;

  -- DESCRIPTION :
  --   Creates an array of doubles with real and imaginary parts of
  --   the complex numbers in the matrix m, stored consecutively 
  --   row after row.

end C_to_Ada_Arrays;
