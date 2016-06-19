with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Dense_Series;             use Standard_Dense_Series;
with Standard_Series_Polynomials;
with Standard_Series_Poly_Systems;

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

  procedure get ( p : out Standard_Series_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure get ( file : in file_type;
                  p : out Standard_Series_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );

  -- DESCRIPTION :
  --   If idx = 0, then all coefficients of p have order zero.
  --   If idx > 0, then the variable of that index is taken to
  --   be the series variable.  If verbose, then extra information
  --   is written to screen during the conversions.

  procedure put ( p : in Standard_Series_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( file : in file_type;
                  p : in Standard_Series_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );

  -- DESCRIPTION :
  --   Converts the series polynomial into a regular polynomial
  --   and writes its representation as an expanded polynomial.
  --   If idx = 0, the order of coefficients of p is taken as zero.
  --   If idx > 0, then the variable with that index is the series variable.
  --   If verbose, then the convertors write extra information to screen.

  procedure get ( ls : out Standard_Series_Poly_Systems.Link_to_Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );

  -- DESCRIPTION :
  --   Reads a polynomial system in several variables, where the variable
  --   with index idx is considered as the series variable if idx > 0.
  --   If idx = 0, then the coefficients of ls on return are series of
  --   order zero.  If verbose, the convertors write information to screen.

  procedure put ( s : in Standard_Series_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( file : in file_type;
                  s : in Standard_Series_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );

  -- DESCRIPTION :
  --   Writes the series polynomials in s as expanded polynomials,
  --   where idx is the index of the series variable, if idx > 0.
  --   If idx = 0, then the coefficients of s have order zero.
  --   If verbose, then the convertors write information to screen.

end Series_and_Polynomials_io;
