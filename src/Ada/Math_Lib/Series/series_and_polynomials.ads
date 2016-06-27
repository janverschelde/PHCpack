with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with Standard_Dense_Series;
with Standard_Dense_Series_Vectors;
with Standard_Dense_Series_Matrices;
with Standard_Series_Polynomials;
with Standard_Series_Poly_Systems;
with DoblDobl_Dense_Series;
with DoblDobl_Dense_Series_Vectors;
with DoblDobl_Dense_Series_Matrices;
with DoblDobl_Series_Polynomials;
with DoblDobl_Series_Poly_Systems;
with QuadDobl_Dense_Series;
with QuadDobl_Dense_Series_Vectors;
with QuadDobl_Dense_Series_Matrices;
with QuadDobl_Series_Polynomials;
with QuadDobl_Series_Poly_Systems;

package Series_and_Polynomials is

-- DESCRIPTION :
--   Exports functions to convert polynomials with complex coefficients
--   into polynomials with series as coefficients, and vice versa.
--   The conversion routines give immediate access to symbolic i/o.
--   Three levels of precision are supported:
--   standard double, double double, and quad double precision.

  function Series_to_Polynomial
             ( s : Standard_Dense_Series.Series )
             return Standard_Complex_Polynomials.Poly;
  function Series_to_Polynomial
             ( s : DoblDobl_Dense_Series.Series )
             return DoblDobl_Complex_Polynomials.Poly;
  function Series_to_Polynomial
             ( s : QuadDobl_Dense_Series.Series )
             return QuadDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns the representation of the series s as a polynomial
  --   in one variable with complex coefficients,
  --   in double, double double, or quad double precision.
  --   This conversion is useful for symbolic output of a series.

  function Polynomial_to_Series
             ( p : Standard_Complex_Polynomials.Poly;
               idx : integer32 := 1 )
             return Standard_Dense_Series.Series;
  function Polynomial_to_Series
             ( p : DoblDobl_Complex_Polynomials.Poly;
               idx : integer32 := 1 )
             return DoblDobl_Dense_Series.Series;
  function Polynomial_to_Series
             ( p : QuadDobl_Complex_Polynomials.Poly;
               idx : integer32 := 1 )
             return QuadDobl_Dense_Series.Series;

  -- DESCRIPTION :
  --   Given in p a polynomial where the series variable has index idx,
  --   with complex coefficients, returns the series representation of p,
  --   in double, double double, or quad double precision.
  --   This conversion is useful for symbolic input of a series.

  -- REQUIRED : degree(p) <= Standard_Dense_Series.max_order.

  function Polynomial_to_Series_Polynomial
             ( p : Standard_Complex_Polynomials.Poly;
               idx : integer32 := 0; verbose : boolean := false )
             return Standard_Series_Polynomials.Poly;
  function Polynomial_to_Series_Polynomial
             ( p : DoblDobl_Complex_Polynomials.Poly;
               idx : integer32 := 0; verbose : boolean := false )
             return DoblDobl_Series_Polynomials.Poly;
  function Polynomial_to_Series_Polynomial
             ( p : QuadDobl_Complex_Polynomials.Poly;
               idx : integer32 := 0; verbose : boolean := false )
             return QuadDobl_Series_Polynomials.Poly;

  -- DESCRIPTION :
  --   By default, if idx is zero, then the coefficient of each term in p
  --   is copied into a series of order zero.  Otherwise, if idx > 0,
  --   then the variable with that index in p is the series parameter.
  --   For example, t^3*x + 2*t*x becomes (2*t + t^3)*x, if idx = 1.
  --   If idx > 0, the number of variables in the polynomial on return is 
  --   one less than the number of variables in the input polynomial p.
  --   Extra output is written to screen if verbose is true.

  function Series_Polynomial_to_Polynomial
             ( s : Standard_Series_Polynomials.Poly;
               idx : integer32 := 0; verbose : boolean := false )
             return Standard_Complex_Polynomials.Poly;
  function Series_Polynomial_to_Polynomial
             ( s : DoblDobl_Series_Polynomials.Poly;
               idx : integer32 := 0; verbose : boolean := false )
             return DoblDobl_Complex_Polynomials.Poly;
  function Series_Polynomial_to_Polynomial
             ( s : QuadDobl_Series_Polynomials.Poly;
               idx : integer32 := 0; verbose : boolean := false )
             return QuadDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Converts a polynomial s with coefficients as series
  --   into a polynomial with complex coefficients.
  --   The first variable in the polynomial on return is
  --   the parameter in the series coefficient of s.
  --   Extra output is written to screen if verbose is true.

  function System_to_Series_System
             ( p : Standard_Complex_Poly_Systems.Poly_Sys;
               idx : integer32 := 0; verbose : boolean := false )
             return Standard_Series_Poly_Systems.Poly_Sys;
  function System_to_Series_System
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               idx : integer32 := 0; verbose : boolean := false )
             return DoblDobl_Series_Poly_Systems.Poly_Sys;
  function System_to_Series_System
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               idx : integer32 := 0; verbose : boolean := false )
             return QuadDobl_Series_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Calls the Polynomial_to_Series_Polynomial to each p(i)
  --   and returns the corresponding system of polynomials
  --   which have as coefficients series with complex coefficients,
  --   in double, double double, or quad double precision.

  function Series_System_to_System
             ( s : Standard_Series_Poly_Systems.Poly_Sys;
               idx : integer32 := 0; verbose : boolean := false )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function Series_System_to_System
             ( s : DoblDobl_Series_Poly_Systems.Poly_Sys;
               idx : integer32 := 0; verbose : boolean := false )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Series_System_to_System
             ( s : QuadDobl_Series_Poly_Systems.Poly_Sys;
               idx : integer32 := 0; verbose : boolean := false )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Calls the Series_Polynomial_to_Polynomial to each p(i)
  --   and returns the corresponding system of series polynomials,
  --   in double, double double, and quad double precision.

  procedure Set_Order ( v : in out Standard_Dense_Series_Vectors.Vector;
                        order : in integer32 );
  procedure Set_Order ( v : in out DoblDobl_Dense_Series_Vectors.Vector;
                        order : in integer32 );
  procedure Set_Order ( v : in out QuadDobl_Dense_Series_Vectors.Vector;
                        order : in integer32 );

  -- DESCRIPTION :
  --   Sets every series in the vector v to the given order,
  --   filling in with zero coefficients if the given order
  --   is larger than the current order.
  --   For evaluation in polynomials, the order must be set high enough
  --   so that the presence of higher degrees is noticed.

  procedure Set_Order ( m : in out Standard_Dense_Series_Matrices.Matrix;
                        order : in integer32 );
  procedure Set_Order ( m : in out DoblDobl_Dense_Series_Matrices.Matrix;
                        order : in integer32 );
  procedure Set_Order ( m : in out QuadDobl_Dense_Series_Matrices.Matrix;
                        order : in integer32 );

  -- DESCRIPTION :
  --   Sets every series in the matrix m to the given order,
  --   filling in with zero coefficients if the given order
  --   is larger than the current order.

  procedure Set_Order ( p : in out Standard_Series_Polynomials.Poly;
                        order : in integer32 );
  procedure Set_Order ( p : in out DoblDobl_Series_Polynomials.Poly;
                        order : in integer32 );
  procedure Set_Order ( p : in out QuadDobl_Series_Polynomials.Poly;
                        order : in integer32 );
  procedure Set_Order ( p : in out Standard_Series_Poly_Systems.Poly_Sys;
                        order : in integer32 );
  procedure Set_Order ( p : in out DoblDobl_Series_Poly_Systems.Poly_Sys;
                        order : in integer32 );
  procedure Set_Order ( p : in out QuadDobl_Series_Poly_Systems.Poly_Sys;
                        order : in integer32 );

  -- DESCRIPTION :
  --   Sets the order of every term in p to the given order.
  --   If the given order is larger than the current order,
  --   then the extra coefficients are set to zero.

  procedure Filter ( s : in out Standard_Dense_Series_Vectors.Vector;
                     tol : in double_float );
  procedure Filter ( s : in out DoblDobl_Dense_Series_Vectors.Vector;
                     tol : in double_float );
  procedure Filter ( s : in out QuadDobl_Dense_Series_Vectors.Vector;
                     tol : in double_float );

  -- DESCRIPTION :
  --   All coefficients in the series of s which are less than tol
  --   in magnitude are set to zero.

end Series_and_Polynomials;
