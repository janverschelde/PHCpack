with text_io;                           use text_io;
with Standard_Dense_Series;             use Standard_Dense_Series;

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

end Series_and_Polynomials_io;
