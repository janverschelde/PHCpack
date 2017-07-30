with text_io;                         use text_io;
with Standard_Integer_Numbers;        use Standard_Integer_Numbers;
with Quad_Double_Numbers;             use Quad_Double_Numbers;
with QuadDobl_Complex_Vectors;        use QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Poly_Systems;   use QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;   use QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Solutions;      use QuadDobl_Complex_Solutions;

package QuadDobl_Hypersurface_Witdrivers is

-- DESCRIPTION :
--   This package provides the drivers to compute a witness set of
--   a hypersurface defined by one polynomial in several variables,
--   in quad double precision.

  function Embedded_System 
               ( n : in integer32;
                 p : in QuadDobl_Complex_Polynomials.Poly;
                 b,v,t : in Vector )
               return Poly_Sys;

  -- DESCRIPTION :
  --   Returns the embedded system in the extrinsic format,
  --   converting the intrinsic line x = b + t*v.

  procedure Call_Root_Finder ( p : in QuadDobl_Complex_Polynomials.Poly );

  -- DESCRIPTION :
  --   Interactive driver to compute a witness set for p
  --   and to write the witness set to file.
  --   The user is prompted for information.

  procedure Call_Root_Finder
               ( file : in file_type;
                 p : in QuadDobl_Complex_Polynomials.Poly;
                 output : in boolean;
                 eps : in quad_double; fail : out boolean );

  -- DESCRIPTION :
  --   Computes a witness set for p and writes it to file.

  -- ON ENTRY :
  --   file      to write the results on;
  --   p         a polynomial in several variables;
  --   output    set to true if intermediate output is desired;
  --   eps       accuracy requirement for univariate root finder.

  -- ON RETURN :
  --   fail      true if accuracy is not met, false otherwise.

  procedure Silent_Root_Finder
               ( p : in QuadDobl_Complex_Polynomials.Poly;
                 eps : in quad_double; fail : out boolean;
                 e : out Link_to_Poly_Sys; esols : out Solution_List );

  -- DESCRIPTION :
  --   Computes a witness set for p, an ordinary polynomial in several
  --   variables, without intermediate output.

  -- ON ENTRY :
  --   p         an ordinary polynomial in several variables;
  --   eps       accuracy requirement for univariate root finder.

  -- ON RETURN :
  --   fail      true if accuracy is not met, false otherwise;
  --   e         embedded polynomial system;
  --   esols     witness points.

  procedure Silent_Root_Finder
               ( p : in QuadDobl_Complex_Laurentials.Poly;
                 eps : in quad_double; fail : out boolean;
                 e : out Link_to_Laur_Sys; esols : out Solution_List );

  -- DESCRIPTION :
  --   Computes a witness set for p, a Laurent polynomial in several
  --   variables, without intermediate output.

  -- ON ENTRY :
  --   p         an ordinary polynomial in several variables;
  --   eps       accuracy requirement for univariate root finder.

  -- ON RETURN :
  --   fail      true if accuracy is not met, false otherwise;
  --   e         embedded polynomial system;
  --   esols     witness points.

end QuadDobl_Hypersurface_Witdrivers;
