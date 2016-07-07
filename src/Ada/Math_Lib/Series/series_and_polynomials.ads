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
with Standard_Dense_Series_VecVecs;
with Standard_Dense_Series_Matrices;
with Standard_Series_Polynomials;
with Standard_Series_Poly_Systems;
with Standard_Series_Jaco_Matrices;
with DoblDobl_Dense_Series;
with DoblDobl_Dense_Series_Vectors;
with DoblDobl_Dense_Series_VecVecs;
with DoblDobl_Dense_Series_Matrices;
with DoblDobl_Series_Polynomials;
with DoblDobl_Series_Poly_Systems;
with DoblDobl_Series_Jaco_Matrices;
with QuadDobl_Dense_Series;
with QuadDobl_Dense_Series_Vectors;
with QuadDobl_Dense_Series_VecVecs;
with QuadDobl_Dense_Series_Matrices;
with QuadDobl_Series_Polynomials;
with QuadDobl_Series_Poly_Systems;
with QuadDobl_Series_Jaco_Matrices;

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

  -- REQUIRED : degree(p) <= Standard_Dense_Series.max_deg.

  function Series_Vector_to_System
             ( v : Standard_Dense_Series_Vectors.Vector )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function Series_Vector_to_System
             ( v : DoblDobl_Dense_Series_Vectors.Vector )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Series_Vector_to_System
             ( v : QuadDobl_Dense_Series_Vectors.Vector )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Returns the representation of all series in v as polynomials
  --   in one variable with complex coefficients,
  --   in double, double double, or quad double precision.
  --   This conversion is useful for symbolic output of a series.

  function System_to_Series_Vector
             ( p : Standard_Complex_Poly_Systems.Poly_Sys;
               idx : integer32 := 1 )
             return Standard_Dense_Series_Vectors.Vector;
  function System_to_Series_Vector
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               idx : integer32 := 1 )
             return DoblDobl_Dense_Series_Vectors.Vector;
  function System_to_Series_Vector
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               idx : integer32 := 1 )
             return QuadDobl_Dense_Series_Vectors.Vector;

  -- DESCRIPTION :
  --   Given in p a system where the series variable has index idx,
  --   with complex coefficients, returns the series representation of p,
  --   in double, double double, or quad double precision.
  --   This conversion is useful for symbolic input of a series.

  -- REQUIRED : degree(p(k)) <= Standard_Dense_Series.max_deg_

  function Series_VecVec_to_System_Array
             ( v : Standard_Dense_Series_VecVecs.VecVec )
             return Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
  function Series_VecVec_to_System_Array
             ( v : DoblDobl_Dense_Series_VecVecs.VecVec )
             return DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
  function Series_VecVec_to_System_Array
             ( v : QuadDobl_Dense_Series_VecVecs.VecVec )
             return QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys;

  -- DESCRIPTION :
  --   Returns the representation of all series vectors in v 
  --   as polynomial systems in one variable with complex coefficients,
  --   in double, double double, or quad double precision.
  --   This conversion is useful for symbolic output of a series.

  function System_Array_to_Series_VecVec
             ( p : Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
               idx : integer32 := 1 )
             return Standard_Dense_Series_VecVecs.VecVec;
  function System_Array_to_Series_VecVec
             ( p : DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
               idx : integer32 := 1 )
             return DoblDobl_Dense_Series_VecVecs.VecVec;
  function System_Array_to_Series_VecVec
             ( p : QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
               idx : integer32 := 1 )
             return QuadDobl_Dense_Series_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Given in p a system array where the series variable has index idx,
  --   with complex coefficients, returns the series representation of p,
  --   in double, double double, or quad double precision.
  --   This conversion is useful for symbolic input of a series.

  -- REQUIRED : degree(p) <= Standard_Dense_Series.max_deg.

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
  --   is copied into a series of degree zero.  Otherwise, if idx > 0,
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

  procedure Set_Degree ( v : in out Standard_Dense_Series_Vectors.Vector;
                         degree : in integer32 );
  procedure Set_Degree ( v : in out DoblDobl_Dense_Series_Vectors.Vector;
                         degree : in integer32 );
  procedure Set_Degree ( v : in out QuadDobl_Dense_Series_Vectors.Vector;
                         degree : in integer32 );

  -- DESCRIPTION :
  --   Sets every series in the vector v to the given degree,
  --   filling in with zero coefficients if the given degree
  --   is larger than the current degree.
  --   For evaluation in polynomials, the degree must be set high enough
  --   so that the presence of higher degrees is noticed.

  procedure Set_Degree ( m : in out Standard_Dense_Series_Matrices.Matrix;
                         degree : in integer32 );
  procedure Set_Degree ( m : in out DoblDobl_Dense_Series_Matrices.Matrix;
                         degree : in integer32 );
  procedure Set_Degree ( m : in out QuadDobl_Dense_Series_Matrices.Matrix;
                         degree : in integer32 );

  -- DESCRIPTION :
  --   Sets every series in the matrix m to the given degree,
  --   filling in with zero coefficients if the given degree
  --   is larger than the current degree.

  procedure Set_Degree ( p : in out Standard_Series_Polynomials.Poly;
                         degree : in integer32 );
  procedure Set_Degree ( p : in out DoblDobl_Series_Polynomials.Poly;
                         degree : in integer32 );
  procedure Set_Degree ( p : in out QuadDobl_Series_Polynomials.Poly;
                         degree : in integer32 );
  procedure Set_Degree ( p : in out Standard_Series_Poly_Systems.Poly_Sys;
                         degree : in integer32 );
  procedure Set_Degree ( p : in out DoblDobl_Series_Poly_Systems.Poly_Sys;
                         degree : in integer32 );
  procedure Set_Degree ( p : in out QuadDobl_Series_Poly_Systems.Poly_Sys;
                         degree : in integer32 );
  procedure Set_Degree ( jm : in out Standard_Series_Jaco_Matrices.Jaco_Mat;
                         degree : in integer32 );
  procedure Set_Degree ( jm : in out DoblDobl_Series_Jaco_Matrices.Jaco_Mat;
                         degree : in integer32 );
  procedure Set_Degree ( jm : in out QuadDobl_Series_Jaco_Matrices.Jaco_Mat;
                         degree : in integer32 );

  -- DESCRIPTION :
  --   Sets the degree of every term in p or jm to the given degree.
  --   If the given degree is larger than the current degree,
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
