with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs;
with Standard_Dense_Series_Vectors;
with Standard_Series_Polynomials;
with Standard_Series_Poly_Systems;
with DoblDobl_Dense_Series_Vectors;
with DoblDobl_Series_Polynomials;
with DoblDobl_Series_Poly_Systems;
with QuadDobl_Dense_Series_Vectors;
with QuadDobl_Series_Polynomials;
with QuadDobl_Series_Poly_Systems;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;

package Series_Polynomial_Gradients is

-- DESCRIPTION :
--   This package exports utilities to define, evaluate, and differentiate
--   polynomials in several variables with power series coefficients.
--   The functions are mainly wrappers and used for testing purposes.

  function Standard_Polynomial
             ( c : Standard_Speelpenning_Convolutions.Convolution_Circuit )
             return Standard_Series_Polynomials.Poly;
  function DoblDobl_Polynomial
             ( c : DoblDobl_Speelpenning_Convolutions.Convolution_Circuit )
             return DoblDobl_Series_Polynomials.Poly;
  function QuadDobl_Polynomial
             ( c : QuadDobl_Speelpenning_Convolutions.Convolution_Circuit )
             return QuadDobl_Series_Polynomials.Poly;

  -- DESCRIPTION :
  --   Makes the polynomial system equivalent to the convolution circuits.

  function Standard_System
             ( c : Standard_Speelpenning_Convolutions.Convolution_Circuits )
             return Standard_Series_Poly_Systems.Poly_Sys;
  function DoblDobl_System
             ( c : DoblDobl_Speelpenning_Convolutions.Convolution_Circuits )
             return DoblDobl_Series_Poly_Systems.Poly_Sys;
  function QuadDobl_System
             ( c : QuadDobl_Speelpenning_Convolutions.Convolution_Circuits )
             return QuadDobl_Series_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Makes the polynomial system equivalent to the convolution circuits.

  function Standard_Product
             ( dim,deg : integer32 )
             return Standard_Series_Polynomials.Poly;
  function DoblDobl_Product
             ( dim,deg : integer32 )
             return DoblDobl_Series_Polynomials.Poly;
  function QuadDobl_Product
             ( dim,deg : integer32 )
             return QuadDobl_Series_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns the product of the first dim variables,
  --   as a polynomial where the coefficients are truncated power series,
  --   truncated to degree deg, with standard double, double double,
  --   or quad double precision coefficients.

  function Standard_Product
             ( deg : integer32;
               xp : Standard_Integer_Vectors.Vector )
             return Standard_Series_Polynomials.Poly;
  function DoblDobl_Product
             ( deg : integer32;
               xp : Standard_Integer_Vectors.Vector )
             return DoblDobl_Series_Polynomials.Poly;
  function QuadDobl_Product
             ( deg : integer32;
               xp : Standard_Integer_Vectors.Vector )
             return QuadDobl_Series_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns the product of the first dim variables, dim = xp'last,
  --   using the exponent vector in xp,
  --   as a polynomial where the coefficients are truncated power series,
  --   truncated to degree deg, with standard double, double double,
  --   or quad double precision coefficients.

  function Standard_Polynomial
             ( dim,deg : integer32;
               xps : Standard_Integer_VecVecs.VecVec )
             return Standard_Series_Polynomials.Poly;
  function Standard_Polynomial
             ( dim : integer32;
               xps : Standard_Integer_VecVecs.VecVec;
               cff : Standard_Dense_Series_Vectors.Vector;
               isxidx : boolean := true )
             return Standard_Series_Polynomials.Poly;
  function DoblDobl_Polynomial
             ( dim,deg : integer32;
               xps : Standard_Integer_VecVecs.VecVec )
             return DoblDobl_Series_Polynomials.Poly;
  function DoblDobl_Polynomial
             ( dim : integer32;
               xps : Standard_Integer_VecVecs.VecVec;
               cff : DoblDobl_Dense_Series_Vectors.Vector;
               isxidx : boolean := true )
             return DoblDobl_Series_Polynomials.Poly;
  function QuadDobl_Polynomial
             ( dim,deg : integer32;
               xps : Standard_Integer_VecVecs.VecVec )
             return QuadDobl_Series_Polynomials.Poly;
  function QuadDobl_Polynomial
             ( dim : integer32;
               xps : Standard_Integer_VecVecs.VecVec;
               cff : QuadDobl_Dense_Series_Vectors.Vector;
               isxidx : boolean := true )
             return QuadDobl_Series_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns the polynomial in dim variables, with exponents in xps,
  --   and optionally, the coefficients in cff,
  --   as a polynomial where the coefficients are truncated power series,
  --   truncated to degree deg, with standard double, double double,
  --   or quad double precision coefficients.
  --   The isxidx indicates if xps is an exponent index vector or not.
  --   If not, then xps holds the actual values of the exponents,
  --   otherwise, if isxidx, then xps holds the indices of the
  --   variables which appear with exponent one in the monomials.

  -- REQUIRED : cff'range = xps'range.

  function Standard_Gradient
             ( p : Standard_Series_Polynomials.Poly;
               x : Standard_Dense_Series_Vectors.Vector )
             return Standard_Dense_Series_Vectors.Vector;
  function DoblDobl_Gradient
             ( p : DoblDobl_Series_Polynomials.Poly;
               x : DoblDobl_Dense_Series_Vectors.Vector )
             return DoblDobl_Dense_Series_Vectors.Vector;
  function QuadDobl_Gradient
             ( p : QuadDobl_Series_Polynomials.Poly;
               x : QuadDobl_Dense_Series_Vectors.Vector )
             return QuadDobl_Dense_Series_Vectors.Vector;

  -- DESCRIPTION :
  --   Evaluates the gradient of p at x, for testing purposes,
  --   in double, double double, or quad double precision.

  function Standard_Series_Coefficients
             ( s : Standard_Dense_Series_Vectors.Vector )
             return Standard_Complex_VecVecs.VecVec;
  function DoblDobl_Series_Coefficients
             ( s : DoblDobl_Dense_Series_Vectors.Vector )
             return DoblDobl_Complex_VecVecs.VecVec;
  function QuadDobl_Series_Coefficients
             ( s : QuadDobl_Dense_Series_Vectors.Vector )
             return QuadDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns the coefficients of the series in the vector of vectors.
  --   The range of the k-th vector is 0..s(k).deg.

end Series_Polynomial_Gradients;
