with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Deca_Double_Numbers;               use Deca_Double_Numbers;
with DecaDobl_Complex_Numbers;          use DecaDobl_Complex_Numbers;
with DecaDobl_Complex_Series;           use DecaDobl_Complex_Series;

package DecaDobl_Complex_Series_Transforms is

-- DESCRIPTION :
--   If the coefficients of a series s(t) grow too large,
--   then with the transform t = c*t, for the proper value of c,
--   the coefficients of s(c*t) become all smaller in modulus than one.
--   The operations, which run in deca double precision, 
--   compute the proper value of the factor c and then
--   transform the coefficient of the series 
--   so all coefficients have modulus one or less.

  procedure Maximum_Coefficient_Modulus
              ( s : in Series;
                idx : out integer32; maxcff : out deca_double );

  -- DESCRIPTION :
  --   Gives a series of degree s.deg > 1, returns the index of the
  --   coefficient with the maximum modulus and this maximum modulus.

  -- ON ENTRY :
  --   s        a truncated power series.

  -- ON RETURN :
  --   idx      the index in 1..s.deg of the coefficient
  --            which has the maximum modulus,
  --            if s.deg is zero, then idx is zero as well;
  --   maxcff   the maximum modulus of the coefficients of the terms
  --            in s with positive power, is -1.0 if s.deg is zero.

  procedure Coefficient_Modulus_Transform
              ( s : in out Series;
                idx : in integer32; maxcff : in deca_double );

  -- DESCRIPTION :
  --   Transforms the coefficients of the series s using the index idx
  --   and modulus maxcff computed by Maximum_Coefficient_Modulus.

  -- ON ENTRY :
  --   s        a truncated power series;
  --   idx      an index in the range 1..s.deg;
  --   maxcff   modulus of the coefficient s.cff(idx).

  -- ON RETURN :
  --   s        coefficients of power k are divided by c^k,
  --            where c = maxcff**(1.0/idx).

  function Scale ( s : DecaDobl_Complex_Series.Series;
                   c : deca_double )
                 return DecaDobl_Complex_Series.Series;
  function Scale ( s : DecaDobl_Complex_Series.Series;
                   c : Complex_Number )
                 return DecaDobl_Complex_Series.Series;

  -- DESCRIPTION :
  --   Returns the series s(c*t).

  procedure Transform ( s : in out Series );

  -- DESCRIPTION :
  --   Transforms the coefficients in s so their modulus is one or less,
  --   using the procedure Coefficient_Modulus_Transform.

end DecaDobl_Complex_Series_Transforms;
