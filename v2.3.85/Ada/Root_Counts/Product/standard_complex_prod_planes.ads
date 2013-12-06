with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Prod_Systems;      use Standard_Complex_Prod_Systems;

package Standard_Complex_Prod_Planes is

-- DESCRIPTION :
--   This package collects some useful routines to create and solve
--   systems of products of linear equations.

  function Create return Prod_Sys;

  -- DESCRIPTION :
  --   Returns a product system from the hyperplane data stored in
  --   the package Random_Product_System.

  -- REQUIRED : Standard_Linear_Product_System.Dimension is nonzero.

  function Linear_Index
             ( d : Standard_Natural_Vectors.Vector ) return integer32;

  -- DESCRIPTION :
  --   Returns 0 if d is zero, otherwise returns k for d(k) = 1
  --   and all other entries are zero.  Returns -1 in all other cases.

  procedure Hyperplane_Coefficients
              ( p : in Poly; fail : out boolean;
                h : out Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Returns in h the coefficients of the hyperplane p,
  --   fail becomes false if p is nonlinear.

  procedure Store ( p : in Prod_Sys; fail : out boolean );

  -- DESCRIPTION :
  --   Stores the factors in p as a linear-product system,
  --   using the internal data in Standard_Linear_Product_System.
  --   Failure is reported if some factor in p is nonlinear.

  function Degrees ( p : Prod_Sys ) return Standard_Natural_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the number of factors in every equation of p.
  --   If every factor is linear, then the vector on return contains
  --   the degree of every equation.

  function Values_at_Hyperplanes
              ( i : natural32; x : Standard_Complex_Vectors.Vector )
              return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 1..n, where n is the number of
  --   hyperplanes of the i-th linear-product polynomial.
  --   The vector on return contains the values of all hyperplanes
  --   evaluated at x.

  function Eval ( i : natural32; x : Standard_Complex_Vectors.Vector )
                return Complex_Number;

  -- DESCRIPTION :
  --   Returns the value of the i-th linear-product polynomial at x.

  function Eval ( x : Standard_Complex_Vectors.Vector )
                return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the value of the linear-product system at x.

  function Gradient ( i : natural32; x : Standard_Complex_Vectors.Vector )
                    return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the value of the gradient of the i-th linear-product
  --   polynomial at x.

  function Jacobian ( x : Standard_Complex_Vectors.Vector ) return Matrix;

  -- DESCRIPTION :
  --   Returns the value of the Jacobian matrix at x.

end Standard_Complex_Prod_Planes;
