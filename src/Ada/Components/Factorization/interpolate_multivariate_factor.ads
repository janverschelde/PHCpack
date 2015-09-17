with text_io;                            use text_io;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;

package Interpolate_Multivariate_Factor is

-- DESCRIPTION :
--   This package provides function to construct a symbolic form of
--   a factor of a polynomial p, given p and witness points on the factor,
--   with interpolation, without and with output to file,
--   in standard double, double double or quad double precision.
--   The polynomials we interpolate are normalized.

  procedure Normalize ( p : in out Standard_Complex_Polynomials.Poly );
  procedure Normalize ( p : in out DoblDobl_Complex_Polynomials.Poly );
  procedure Normalize ( p : in out QuadDobl_Complex_Polynomials.Poly );
  procedure Normalize ( p : in out Standard_Complex_Poly_Systems.Poly_Sys );
  procedure Normalize ( p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Normalize ( p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys );

  -- DESCRIPTION :
  --   A polynomial is normalized if its leading coefficient is one.
  --   This procedure normalizes the polynomial dividing the polynomial
  --   by its leading coefficient.

  function Interpolate_Factor
              ( p : Standard_Complex_Polynomials.Poly;
                b,v,w : Standard_Complex_Vectors.Vector )
              return Standard_Complex_Polynomials.Poly;
  function Interpolate_Factor
              ( p : DoblDobl_Complex_Polynomials.Poly;
                b,v,w : DoblDobl_Complex_Vectors.Vector )
              return DoblDobl_Complex_Polynomials.Poly;
  function Interpolate_Factor
              ( p : QuadDobl_Complex_Polynomials.Poly;
                b,v,w : QuadDobl_Complex_Vectors.Vector )
              return QuadDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Finds a polynomial interpolating through a factor of p,
  --   without any intermediate output.

  -- NOTICE : 
  --   For n = 2, the Newton-Taylor form is used, while for n > 2,
  --   we rely on traces to find the expanded form.

  function Interpolate_Factor
              ( file : file_type;
                p : Standard_Complex_Polynomials.Poly;
                b,v,w : Standard_Complex_Vectors.Vector )
              return Standard_Complex_Polynomials.Poly;
  function Interpolate_Factor
              ( file : file_type;
                p : DoblDobl_Complex_Polynomials.Poly;
                b,v,w : DoblDobl_Complex_Vectors.Vector )
              return DoblDobl_Complex_Polynomials.Poly;
  function Interpolate_Factor
              ( file : file_type;
                p : QuadDobl_Complex_Polynomials.Poly;
                b,v,w : QuadDobl_Complex_Vectors.Vector )
              return QuadDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Finds a polynomial interpolating through a factor of p,
  --   with intermediate output.

end Interpolate_Multivariate_Factor;
