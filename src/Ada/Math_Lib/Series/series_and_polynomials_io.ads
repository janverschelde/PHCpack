with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Dense_Series;
with Standard_Dense_Series_Vectors;
with Standard_Series_Polynomials;
with Standard_Series_Poly_Systems;
with DoblDobl_Dense_Series;
with DoblDobl_Dense_Series_Vectors;
with DoblDobl_Series_Polynomials;
with DoblDobl_Series_Poly_Systems;
with QuadDobl_Dense_Series;
with QuadDobl_Dense_Series_Vectors;
with QuadDobl_Series_Polynomials;
with QuadDobl_Series_Poly_Systems;

package Series_and_Polynomials_io is

-- DESCRIPTION :
--   Provides symbolic input and output of truncated power series
--   and polynomials with series coefficients.

  procedure get ( s : out Standard_Dense_Series.Series );
  procedure get ( s : out DoblDobl_Dense_Series.Series );
  procedure get ( s : out QuadDobl_Dense_Series.Series );
  procedure get ( file : in file_type;
                  s : out Standard_Dense_Series.Series );
  procedure get ( file : in file_type;
                  s : out DoblDobl_Dense_Series.Series );
  procedure get ( file : in file_type;
                  s : out QuadDobl_Dense_Series.Series );

  -- DESCRIPTION :
  --   Reads a polynomial in one variable and returns the
  --   corresponding series representation.

  procedure put ( s : in Standard_Dense_Series.Series );
  procedure put ( s : in DoblDobl_Dense_Series.Series );
  procedure put ( s : in QuadDobl_Dense_Series.Series );
  procedure put ( file : in file_type;
                  s : in Standard_Dense_Series.Series );
  procedure put ( file : in file_type;
                  s : in DoblDobl_Dense_Series.Series );
  procedure put ( file : in file_type;
                  s : in QuadDobl_Dense_Series.Series );

  -- DESCRIPTION :
  --   Writes a polynomial in one variable,
  --   representing the series in s.

  procedure get ( lv : out Standard_Dense_Series_Vectors.Link_to_Vector;
                  idx : in integer32 := 1; verbose : in boolean := false );
  procedure get ( lv : out DoblDobl_Dense_Series_Vectors.Link_to_Vector;
                  idx : in integer32 := 1; verbose : in boolean := false );
  procedure get ( lv : out QuadDobl_Dense_Series_Vectors.Link_to_Vector;
                  idx : in integer32 := 1; verbose : in boolean := false );

  -- DESCRIPTION :
  --   Asks the user first if the series are on file.
  --   Then reads a number of univariate polynomials and converts
  --   this sequence of polynomials into a vector of series.
  --   The format on file is the same as a polynomial system.
  --   By default, the series variable is assumed to be the first variable
  --   of the polynomials read.  Given a value to idx changes this default.

  procedure put ( v : in Standard_Dense_Series_Vectors.Vector );
  procedure put ( v : in DoblDobl_Dense_Series_Vectors.Vector );
  procedure put ( v : in QuadDobl_Dense_Series_Vectors.Vector );
  procedure put ( file : in file_type;
                  v : in Standard_Dense_Series_Vectors.Vector );
  procedure put ( file : in file_type;
                  v : in DoblDobl_Dense_Series_Vectors.Vector );
  procedure put ( file : in file_type;
                  v : in QuadDobl_Dense_Series_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes the series as a system of v'length univariate polynomials
  --   to standard output or to file.  The first line contains v'length,
  --   followed by one, in the same formate as accepted by get.

  procedure get ( p : out Standard_Series_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure get ( p : out DoblDobl_Series_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure get ( p : out QuadDobl_Series_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure get ( file : in file_type;
                  p : out Standard_Series_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure get ( file : in file_type;
                  p : out DoblDobl_Series_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure get ( file : in file_type;
                  p : out QuadDobl_Series_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );

  -- DESCRIPTION :
  --   If idx = 0, then all coefficients of p have degree zero.
  --   If idx > 0, then the variable of that index is taken to
  --   be the series variable.  If verbose, then extra information
  --   is written to screen during the conversions.

  procedure put ( p : in Standard_Series_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( p : in DoblDobl_Series_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( p : in QuadDobl_Series_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( file : in file_type;
                  p : in Standard_Series_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( file : in file_type;
                  p : in DoblDobl_Series_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( file : in file_type;
                  p : in QuadDobl_Series_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );

  -- DESCRIPTION :
  --   Converts the series polynomial into a regular polynomial
  --   and writes its representation as an expanded polynomial.
  --   If idx = 0, the degree of coefficients of p is taken as zero.
  --   If idx > 0, then the variable with that index is the series variable.
  --   If verbose, then the convertors write extra information to screen.

  procedure get ( ls : out Standard_Series_Poly_Systems.Link_to_Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure get ( ls : out DoblDobl_Series_Poly_Systems.Link_to_Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure get ( ls : out QuadDobl_Series_Poly_Systems.Link_to_Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );

  -- DESCRIPTION :
  --   First asks the user whether the system is on file.
  --   Reads a polynomial system in several variables, where the variable
  --   with index idx is considered as the series variable if idx > 0.
  --   If idx = 0, then the coefficients of ls on return are series of
  --   degree zero.  If verbose, the convertors write information to screen.

  procedure put ( s : in Standard_Series_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( s : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( s : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( file : in file_type;
                  s : in Standard_Series_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( file : in file_type;
                  s : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( file : in file_type;
                  s : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );

  -- DESCRIPTION :
  --   Writes the series polynomials in s as expanded polynomials,
  --   where idx is the index of the series variable, if idx > 0.
  --   If idx = 0, then the coefficients of s have degree zero.
  --   If verbose, then the convertors write information to screen.

end Series_and_Polynomials_io;
