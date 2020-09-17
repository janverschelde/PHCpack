with Deca_Double_Numbers;               use Deca_Double_Numbers;
with DecaDobl_Complex_Numbers;          use DecaDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with DecaDobl_Complex_Vectors;
with DecaDobl_Complex_Series_Vectors;
with DecaDobl_Complex_Series_VecVecs;

package DecaDobl_CSeries_Vector_Functions is

-- DESCRIPTION :
--   Functions to evaluate vectors of power series,
--   in deca double precision.

  function Eval ( v : DecaDobl_Complex_Series_Vectors.Vector;
                  t : deca_double )
                return DecaDobl_Complex_Vectors.Vector;
  function Eval ( v : DecaDobl_Complex_Series_Vectors.Vector;
                  t : Complex_Number )
                return DecaDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the values of all series in v at t.

  function Eval ( v : DecaDobl_Complex_Series_Vectors.Vector;
                  w : Standard_Integer_Vectors.Vector;
                  t : deca_double )
                return DecaDobl_Complex_Vectors.Vector;
  function Eval ( v : DecaDobl_Complex_Series_Vectors.Vector;
                  w : Standard_Integer_Vectors.Vector;
                  t : Complex_Number )
                return DecaDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the values of all series in v at t,
  --   with weighted powers w(k)/w(w'last) for the k-th series.

  -- REQUIRED : w'range = v'first..v'last+1
  --   and t /= 0 if there are negative weights in w.

  function Shift ( v : DecaDobl_Complex_Series_Vectors.Vector;
                   c : deca_double )
                 return DecaDobl_Complex_Series_Vectors.Vector;
  function Shift ( v : DecaDobl_Complex_Series_Vectors.Vector;
                   c : Complex_Number )
                 return DecaDobl_Complex_Series_Vectors.Vector;
  function Shift ( v : DecaDobl_Complex_Series_VecVecs.VecVec;
                   c : deca_double )
                 return DecaDobl_Complex_Series_VecVecs.VecVec;
  function Shift ( v : DecaDobl_Complex_Series_VecVecs.VecVec;
                   c : Complex_Number )
                 return DecaDobl_Complex_Series_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns the series v, shifted with the constant c.
  --   The returned series coefficients place the center of the
  --   evaluation from c to 0, such that evaluating the returned
  --   series at 0 is the same as replacing the series parameter t by c.

  procedure Shift ( v : in out DecaDobl_Complex_Series_Vectors.Vector;
                    c : in deca_double );
  procedure Shift ( v : in out DecaDobl_Complex_Series_Vectors.Vector;
                    c : in Complex_Number );
  procedure Shift ( v : in out DecaDobl_Complex_Series_VecVecs.VecVec;
                    c : in deca_double );
  procedure Shift ( v : in out DecaDobl_Complex_Series_VecVecs.VecVec;
                    c : in Complex_Number );

  -- DESCRIPTION :
  --   Shifts all series in v with the constant c, changing the series
  --   coefficients with the same effect that the series parameter t
  --   is replaced by t-c, so that Eval(v,c) = Eval(Shift(v,c),0).

  function Make_Deep_Copy
             ( v : DecaDobl_Complex_Series_Vectors.Vector )
             return DecaDobl_Complex_Series_Vectors.Vector;
  function Make_Deep_Copy
             ( v : DecaDobl_Complex_Series_VecVecs.VecVec )
             return DecaDobl_Complex_Series_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns a deep copy of all series in v.

  procedure Deep_Clear ( v : in out DecaDobl_Complex_Series_Vectors.Vector );
  procedure Deep_Clear ( v : in out DecaDobl_Complex_Series_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Deallocates all memory occupied by the series in v.

end DecaDobl_CSeries_Vector_Functions;
