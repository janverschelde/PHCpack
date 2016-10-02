with Double_Double_Numbers;             use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;          use DoblDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with DoblDobl_Complex_Vectors;
with DoblDobl_Dense_Series_Vectors;

package DoblDobl_Series_Vector_Functions is

-- DESCRIPTION :
--   Functions to evaluate vectors of power series.

  function Eval ( v : DoblDobl_Dense_Series_Vectors.Vector;
                  t : double_double )
                return DoblDobl_Complex_Vectors.Vector;
  function Eval ( v : DoblDobl_Dense_Series_Vectors.Vector;
                  t : Complex_Number )
                return DoblDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the values of all series in v at t.

  function Eval ( v : DoblDobl_Dense_Series_Vectors.Vector;
                  w : Standard_Integer_Vectors.Vector;
                  t : double_double )
                return DoblDobl_Complex_Vectors.Vector;
  function Eval ( v : DoblDobl_Dense_Series_Vectors.Vector;
                  w : Standard_Integer_Vectors.Vector;
                  t : Complex_Number )
                return DoblDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the values of all series in v at t,
  --   with weighted powers w(k)/w(w'last) for the k-th series.

  -- REQUIRED : w'range = v'first..v'last+1
  --   and t /= 0 if there are negative weights in w.

end DoblDobl_Series_Vector_Functions;
