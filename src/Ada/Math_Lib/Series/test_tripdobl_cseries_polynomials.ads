with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with TripDobl_Complex_Series;
with TripDobl_Complex_Series_Vectors;
with TripDobl_CSeries_Polynomials;

package Test_TripDobl_CSeries_Polynomials is

-- DESCRIPTION :
--   Tests polynomials in several variables, with as coefficients
--   truncated power series in triple double precision.
 
  procedure Write ( s : in TripDobl_CSeries_Polynomials.Poly );

  -- DESCRIPTION :
  --   Very simple output of a polynomial with series coefficients.
 
  procedure TripDobl_Test_Conversion;

  -- DESCRIPTION :
  --   Prompts for the number of variables.
  --   and reads in a regular polynomial in several variables,
  --   for conversion into a series polynomial.

  function Factor ( n,k : integer32;
                    s : TripDobl_Complex_Series.Series )
                  return TripDobl_CSeries_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns x[k] - s as a polynomial in n variables.
  --   All coefficients of the polynomial on return have the same degree,
  --   the same as s.deg.

  -- REQUIRED : k is in range 1..n.

  function Product ( s : TripDobl_Complex_Series_Vectors.Vector )
                   return TripDobl_CSeries_Polynomials.Poly;

  -- DESCRIPION :
  --   Returns the product of the factors x[k] - s[k], for k in s'range,
  --   where s'first = 1, in triple double precision.

  procedure TripDobl_Frequent_Evaluation 
              ( dim,deg : in integer32;
                pol : in TripDobl_CSeries_Polynomials.Poly );

  -- DESCRIPTION :
  --   Prompts the user for a number and the evaluates the
  --   polynomial as many times as the number at a random series.
  --   Memory consumption during the loop should remain constant.

  -- ON ENTRY :
  --   dim      number of variables in the polynomial pol;
  --   deg      degree of the power series;
  --   pol      a polynomial in dim variables with power series
  --            coefficients of degree deg.

  procedure TripDobl_Test_Evaluation;

  -- DESCRIPTION :
  --   Prompts for the number of variables and the degree of the series.
  --   Then as many random series as the number of variables are generated.
  --   The polynomial is of the product of x[k] - s[k], where k ranges
  --   over the number of variables 'x' and series 's'.
  --   So the evaluation at the series should produce zero.

  procedure TripDobl_Test_Input_Output;

  -- DESCRIPTION :
  --   Reads a series in symbolic format and writes the series back.

  procedure TripDobl_Test_Polynomial_Series;

  -- DESCRIPTION :
  --   Reads a polynomial in several variables and converts
  --   the polynomial into a polynomial series.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_TripDobl_CSeries_Polynomials;
