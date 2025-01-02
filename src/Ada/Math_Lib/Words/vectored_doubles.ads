with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Vectors;

package Vectored_Doubles is

-- DESCRIPTION :
--   The vectored product of quartered doubles provides a baseline
--   to compare the cost overhead of double double arithmetic 
--   to the cost of double arithmetic.

  procedure Balanced_Quarter_Product
              ( dim : in integer32;
                x0,x1,x2,x3 : in Standard_Floating_Vectors.Vector;
                y0,y1,y2,y3 : in Standard_Floating_Vectors.Vector;
                s0,s1,s2,s3 : out double_float );

  -- DESCRIPTION :
  --   Given balanced quarter doubles in vectors of length dim,
  --   returns the subsums of their inner product.

  -- ON ENTRY :
  --   dim      all vectors have range 1..dim;
  --   x0, ..., x3 are the balanced quarter doubles of x;
  --   y0, ..., y3 are the balanced quarter doubles of y.

  -- ON RETURN :
  --   s0, ..., s3 are the subsums of the inner product of x and y.

  procedure Write_Subsums ( s0,s1,s2,s3 : in double_float );

  -- DESCRIPTION :
  --   Writes all subsums of the convolutions of a product of two vectors.

  function to_double ( s0,s1,s2,s3 : double_float;
                       verbose : boolean := true ) return double_float;

  -- DESCRIPTION :
  --   Adds the subsums in s0, s1, s2, s3 to the return double.
  --   If verbose, then the subsums are written.

end Vectored_Doubles; 
