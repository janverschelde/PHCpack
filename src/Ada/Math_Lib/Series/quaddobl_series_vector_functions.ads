with Quad_Double_Numbers;               use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;          use QuadDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with QuadDobl_Complex_Vectors;
with QuadDobl_Dense_Series_Vectors;

package QuadDobl_Series_Vector_Functions is

-- DESCRIPTION :
--   Functions to evaluate vectors of power series.

  function Eval ( v : QuadDobl_Dense_Series_Vectors.Vector;
                  t : quad_double )
                return QuadDobl_Complex_Vectors.Vector;
  function Eval ( v : QuadDobl_Dense_Series_Vectors.Vector;
                  t : Complex_Number )
                return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the values of all series in v at t.

  function Eval ( v : QuadDobl_Dense_Series_Vectors.Vector;
                  w : Standard_Integer_Vectors.Vector;
                  t : quad_double )
                return QuadDobl_Complex_Vectors.Vector;
  function Eval ( v : QuadDobl_Dense_Series_Vectors.Vector;
                  w : Standard_Integer_Vectors.Vector;
                  t : Complex_Number )
                return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the values of all series in v at t,
  --   with weighted powers w(k)/w(w'last) for the k-th series.

  -- REQUIRED : w'range = v'first..v'last+1
  --   and t /= 0 if there are negative weights in w.

end QuadDobl_Series_Vector_Functions;
