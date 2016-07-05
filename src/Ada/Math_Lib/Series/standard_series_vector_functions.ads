with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
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

end Standard_Series_Vector_Functions;
