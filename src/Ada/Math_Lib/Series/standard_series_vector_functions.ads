with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Dense_Series_Vectors;

package Standard_Series_Vector_Functions is

-- DESCRIPTION :
--   Functions to evaluate vectors of power series.

  function Eval ( v : Standard_Dense_Series_Vectors.Vector;
                  t : double_float )
                return Standard_Complex_Vectors.Vector;
  function Eval ( v : Standard_Dense_Series_Vectors.Vector;
                  t : Complex_Number )
                return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the values of all series in v at t.

  function Eval ( v : Standard_Dense_Series_Vectors.Vector;
                  w : Standard_Integer_Vectors.Vector;
                  t : double_float )
                return Standard_Complex_Vectors.Vector;
  function Eval ( v : Standard_Dense_Series_Vectors.Vector;
                  w : Standard_Integer_Vectors.Vector;
                  t : Complex_Number )
                return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the values of all series in v at t,
  --   with weighted powers w(k)/w(w'last) for the k-th series.

  -- REQUIRED : w'range = v'first..v'last+1
  --   and t /= 0 if there are negative weights in w.

end Standard_Series_Vector_Functions;
