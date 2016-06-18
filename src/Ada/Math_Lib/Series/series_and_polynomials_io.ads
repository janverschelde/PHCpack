with text_io;                           use text_io;
with Standard_Dense_Series;             use Standard_Dense_Series;
with Standard_Series_Polynomials;

package Series_and_Polynomials_io is

-- DESCRIPTION :
--   Provides symbolic input and output of truncated power series
--   and polynomials with series coefficients.

  procedure get ( s : out Series );
  procedure get ( file : in file_type; s : out Series );

  -- DESCRIPTION :
  --   Reads a polynomial in one variable and returns the
  --   corresponding series representation.

  procedure put ( s : in Series );
  procedure put ( file : in file_type; s : in Series );

  -- DESCRIPTION :
  --   Writes a polynomial in one variable,
  --   representing the series in s.

  procedure get ( p : out Standard_Series_Polynomials.Poly );
  procedure get ( file : in file_type;
                  p : out Standard_Series_Polynomials.Poly );

  -- DESCRIPTION :
  --   Reads a polynomial in n+1 variables and converts this
  --   multivariate polynomial into a series polynomial.

  procedure put ( p : in Standard_Series_Polynomials.Poly );
  procedure put ( file : in file_type;
                  p : in Standard_Series_Polynomials.Poly );

  -- DESCRIPTION :
  --   Converts the series polynomial into a regular polynomial
  --   and writes its representation as an expanded polynomial.

end Series_and_Polynomials_io;
