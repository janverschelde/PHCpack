with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Multprec_Complex_Vectors;          use Multprec_Complex_Vectors;
with Multprec_Complex_Solutions;        use Multprec_Complex_Solutions;

package Extrapolate_Solution_Clusters is

  procedure Extrapolate
              ( file : in file_type;
                order,m,n : in integer32; size : in natural32;
                sols : in Solution_Array; extra : out Vector );

  -- DECRIPTION :
  --   Performs extrapolation on the array of solutions.

  -- ON ENTRY :
  --   file     to write intermediate output and diagnostics;
  --   order    order of the extrapolator;
  --   m        multiplicity;
  --   n        dimension of solution lists;
  --   size     size of the numbers;
  --   sols     is of range 0..order and contains approximations 
  --            for increasing values of the continuation parameter t.

  -- ON RETURN :
  --   extra    extrapolated solution vector.

end Extrapolate_Solution_Clusters;
