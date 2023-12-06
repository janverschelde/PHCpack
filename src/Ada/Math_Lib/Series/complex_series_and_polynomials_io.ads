with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Series;
with Standard_Complex_Series_Vectors;
with Standard_CSeries_Polynomials;
with Standard_CSeries_Poly_Systems;
with DoblDobl_Complex_Series;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_CSeries_Polynomials;
with DoblDobl_CSeries_Poly_Systems;
with TripDobl_Complex_Series;
with TripDobl_Complex_Series_Vectors;
with TripDobl_CSeries_Polynomials;
with TripDobl_CSeries_Poly_Systems;
with QuadDobl_Complex_Series;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_CSeries_Polynomials;
with QuadDobl_CSeries_Poly_Systems;
with PentDobl_Complex_Series;
with PentDobl_Complex_Series_Vectors;
with PentDobl_CSeries_Polynomials;
with PentDobl_CSeries_Poly_Systems;
with OctoDobl_Complex_Series;
with OctoDobl_Complex_Series_Vectors;
with OctoDobl_CSeries_Polynomials;
with OctoDobl_CSeries_Poly_Systems;
with DecaDobl_Complex_Series;
with DecaDobl_Complex_Series_Vectors;
with DecaDobl_CSeries_Polynomials;
with DecaDobl_CSeries_Poly_Systems;
with HexaDobl_Complex_Series;
with HexaDobl_Complex_Series_Vectors;
with HexaDobl_CSeries_Polynomials;
with HexaDobl_CSeries_Poly_Systems;

package Complex_Series_and_Polynomials_io is

-- DESCRIPTION :
--   Provides symbolic input and output of truncated power series
--   and polynomials with series coefficients.

  procedure get ( s : out Standard_Complex_Series.Series );
  procedure get ( s : out DoblDobl_Complex_Series.Series );
  procedure get ( s : out TripDobl_Complex_Series.Series );
  procedure get ( s : out QuadDobl_Complex_Series.Series );
  procedure get ( s : out PentDobl_Complex_Series.Series );
  procedure get ( s : out OctoDobl_Complex_Series.Series );
  procedure get ( s : out DecaDobl_Complex_Series.Series );
  procedure get ( s : out HexaDobl_Complex_Series.Series );
  procedure get ( file : in file_type;
                  s : out Standard_Complex_Series.Series );
  procedure get ( file : in file_type;
                  s : out DoblDobl_Complex_Series.Series );
  procedure get ( file : in file_type;
                  s : out TripDobl_Complex_Series.Series );
  procedure get ( file : in file_type;
                  s : out QuadDobl_Complex_Series.Series );
  procedure get ( file : in file_type;
                  s : out PentDobl_Complex_Series.Series );
  procedure get ( file : in file_type;
                  s : out OctoDobl_Complex_Series.Series );
  procedure get ( file : in file_type;
                  s : out DecaDobl_Complex_Series.Series );
  procedure get ( file : in file_type;
                  s : out HexaDobl_Complex_Series.Series );

  -- DESCRIPTION :
  --   Reads a polynomial in one variable and returns the
  --   corresponding series representation.

  procedure put ( s : in Standard_Complex_Series.Series );
  procedure put ( s : in DoblDobl_Complex_Series.Series );
  procedure put ( s : in TripDobl_Complex_Series.Series );
  procedure put ( s : in QuadDobl_Complex_Series.Series );
  procedure put ( s : in PentDobl_Complex_Series.Series );
  procedure put ( s : in OctoDobl_Complex_Series.Series );
  procedure put ( s : in DecaDobl_Complex_Series.Series );
  procedure put ( s : in HexaDobl_Complex_Series.Series );
  procedure put ( file : in file_type;
                  s : in Standard_Complex_Series.Series );
  procedure put ( file : in file_type;
                  s : in DoblDobl_Complex_Series.Series );
  procedure put ( file : in file_type;
                  s : in TripDobl_Complex_Series.Series );
  procedure put ( file : in file_type;
                  s : in QuadDobl_Complex_Series.Series );
  procedure put ( file : in file_type;
                  s : in PentDobl_Complex_Series.Series );
  procedure put ( file : in file_type;
                  s : in OctoDobl_Complex_Series.Series );
  procedure put ( file : in file_type;
                  s : in DecaDobl_Complex_Series.Series );
  procedure put ( file : in file_type;
                  s : in HexaDobl_Complex_Series.Series );

  -- DESCRIPTION :
  --   Writes a polynomial in one variable,
  --   representing the series in s.

  procedure get ( lv : out Standard_Complex_Series_Vectors.Link_to_Vector;
                  idx : in integer32 := 1; verbose : in boolean := false );
  procedure get ( lv : out DoblDobl_Complex_Series_Vectors.Link_to_Vector;
                  idx : in integer32 := 1; verbose : in boolean := false );
  procedure get ( lv : out TripDobl_Complex_Series_Vectors.Link_to_Vector;
                  idx : in integer32 := 1; verbose : in boolean := false );
  procedure get ( lv : out QuadDobl_Complex_Series_Vectors.Link_to_Vector;
                  idx : in integer32 := 1; verbose : in boolean := false );
  procedure get ( lv : out PentDobl_Complex_Series_Vectors.Link_to_Vector;
                  idx : in integer32 := 1; verbose : in boolean := false );
  procedure get ( lv : out OctoDobl_Complex_Series_Vectors.Link_to_Vector;
                  idx : in integer32 := 1; verbose : in boolean := false );
  procedure get ( lv : out DecaDobl_Complex_Series_Vectors.Link_to_Vector;
                  idx : in integer32 := 1; verbose : in boolean := false );
  procedure get ( lv : out HexaDobl_Complex_Series_Vectors.Link_to_Vector;
                  idx : in integer32 := 1; verbose : in boolean := false );

  -- DESCRIPTION :
  --   Asks the user first if the series are on file.
  --   Then reads a number of univariate polynomials and converts
  --   this sequence of polynomials into a vector of series.
  --   The format on file is the same as a polynomial system.
  --   By default, the series variable is assumed to be the first variable
  --   of the polynomials read.  Given a value to idx changes this default.

  procedure put ( v : in Standard_Complex_Series_Vectors.Vector );
  procedure put ( v : in DoblDobl_Complex_Series_Vectors.Vector );
  procedure put ( v : in TripDobl_Complex_Series_Vectors.Vector );
  procedure put ( v : in QuadDobl_Complex_Series_Vectors.Vector );
  procedure put ( v : in PentDobl_Complex_Series_Vectors.Vector );
  procedure put ( v : in OctoDobl_Complex_Series_Vectors.Vector );
  procedure put ( v : in DecaDobl_Complex_Series_Vectors.Vector );
  procedure put ( v : in HexaDobl_Complex_Series_Vectors.Vector );
  procedure put ( file : in file_type;
                  v : in Standard_Complex_Series_Vectors.Vector );
  procedure put ( file : in file_type;
                  v : in DoblDobl_Complex_Series_Vectors.Vector );
  procedure put ( file : in file_type;
                  v : in TripDobl_Complex_Series_Vectors.Vector );
  procedure put ( file : in file_type;
                  v : in QuadDobl_Complex_Series_Vectors.Vector );
  procedure put ( file : in file_type;
                  v : in PentDobl_Complex_Series_Vectors.Vector );
  procedure put ( file : in file_type;
                  v : in OctoDobl_Complex_Series_Vectors.Vector );
  procedure put ( file : in file_type;
                  v : in DecaDobl_Complex_Series_Vectors.Vector );
  procedure put ( file : in file_type;
                  v : in HexaDobl_Complex_Series_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes the series as a system of v'length univariate polynomials
  --   to standard output or to file.  The first line contains v'length,
  --   followed by one, in the same formate as accepted by get.

  procedure get ( p : out Standard_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure get ( p : out DoblDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure get ( p : out TripDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure get ( p : out QuadDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure get ( p : out PentDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure get ( p : out OctoDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure get ( p : out DecaDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure get ( p : out HexaDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure get ( file : in file_type;
                  p : out Standard_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure get ( file : in file_type;
                  p : out DoblDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure get ( file : in file_type;
                  p : out TripDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure get ( file : in file_type;
                  p : out QuadDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure get ( file : in file_type;
                  p : out PentDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure get ( file : in file_type;
                  p : out OctoDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure get ( file : in file_type;
                  p : out DecaDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure get ( file : in file_type;
                  p : out HexaDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );

  -- DESCRIPTION :
  --   If idx = 0, then all coefficients of p have degree zero.
  --   If idx > 0, then the variable of that index is taken to
  --   be the series variable.  If verbose, then extra information
  --   is written to screen during the conversions.

  procedure put ( p : in Standard_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( p : in DoblDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( p : in TripDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( p : in QuadDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( p : in PentDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( p : in OctoDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( p : in DecaDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( p : in HexaDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( file : in file_type;
                  p : in Standard_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( file : in file_type;
                  p : in DoblDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( file : in file_type;
                  p : in TripDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( file : in file_type;
                  p : in QuadDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( file : in file_type;
                  p : in PentDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( file : in file_type;
                  p : in OctoDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( file : in file_type;
                  p : in DecaDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( file : in file_type;
                  p : in HexaDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false );

  -- DESCRIPTION :
  --   Converts the series polynomial into a regular polynomial
  --   and writes its representation as an expanded polynomial.
  --   If idx = 0, the degree of coefficients of p is taken as zero.
  --   If idx > 0, then the variable with that index is the series variable.
  --   If verbose, then the convertors write extra information to screen.

  procedure get ( ls : out Standard_CSeries_Poly_Systems.Link_to_Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure get ( ls : out DoblDobl_CSeries_Poly_Systems.Link_to_Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure get ( ls : out TripDobl_CSeries_Poly_Systems.Link_to_Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure get ( ls : out QuadDobl_CSeries_Poly_Systems.Link_to_Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure get ( ls : out PentDobl_CSeries_Poly_Systems.Link_to_Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure get ( ls : out OctoDobl_CSeries_Poly_Systems.Link_to_Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure get ( ls : out DecaDobl_CSeries_Poly_Systems.Link_to_Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure get ( ls : out HexaDobl_CSeries_Poly_Systems.Link_to_Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );

  -- DESCRIPTION :
  --   First asks the user whether the system is on file.
  --   Reads a polynomial system in several variables, where the variable
  --   with index idx is considered as the series variable if idx > 0.
  --   If idx = 0, then the coefficients of ls on return are series of
  --   degree zero.  If verbose, the convertors write information to screen.

  procedure put ( s : in Standard_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( s : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( s : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( s : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( s : in PentDobl_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( s : in OctoDobl_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( s : in DecaDobl_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( s : in HexaDobl_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( file : in file_type;
                  s : in Standard_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( file : in file_type;
                  s : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( file : in file_type;
                  s : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( file : in file_type;
                  s : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( file : in file_type;
                  s : in PentDobl_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( file : in file_type;
                  s : in OctoDobl_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( file : in file_type;
                  s : in DecaDobl_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );
  procedure put ( file : in file_type;
                  s : in HexaDobl_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false );

  -- DESCRIPTION :
  --   Writes the series polynomials in s as expanded polynomials,
  --   where idx is the index of the series variable, if idx > 0.
  --   If idx = 0, then the coefficients of s have degree zero.
  --   If verbose, then the convertors write information to screen.

end Complex_Series_and_Polynomials_io;
