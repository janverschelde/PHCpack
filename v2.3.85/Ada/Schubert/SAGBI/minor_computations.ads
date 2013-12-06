with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Floating_Matrices;         use Standard_Floating_Matrices;

package Minor_Computations is

  procedure Minors ( file : in file_type;
                     n,m : in natural32; mat : in Matrix );

  -- DESCRIPTION :
  --   Computes all d-minors of a (nxm)-matrix mat, m < n, for d=2,..,m
  --   and writes the results on file.

  function Number_of_Minors ( n,m : natural32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of all d-minors of an n-by-m matrix, for d=2,..,m.

  function Sign_of_Minors ( n,m : natural32; mat : Matrix ) return Vector;

  -- DESCRIPTION :
  --   Returns the sign pattern of all minors in the vector on return.
  --   Every entry of that vector contains either -1, 0, or +1.

end Minor_Computations;
