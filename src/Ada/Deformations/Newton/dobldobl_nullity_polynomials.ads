with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Vectors;          use Standard_Natural_Vectors;
with DoblDobl_Complex_Vectors;          use DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Polynomials;      use DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;     use DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Matrices;
 
package DoblDobl_Nullity_Polynomials is

-- DESCRIPTION :
--   This package collects the auxiliary operations necessary to build
--   and evaluate nullity matrices using double double complex arithmetic.
--   Polynomials subjected to the operations are "nullity polynomials".

  function Derivative 
              ( p : Poly; m : Standard_Natural_Vectors.Vector ) return Poly;

  -- DESCRIPTION :
  --   Returns the derivative of p with respect to the variables in m.

  function Factorial ( m : Standard_Natural_Vectors.Vector ) return natural32;

  -- DESCRIPTION :
  --   Returns the product of all factorials of the elements in m.

  function Monomial_Multiple
              ( m : Standard_Natural_Vectors.Vector;
                f : Poly_Sys ) return Poly_Sys;

  -- DESCRIPTION :
  --   Returns x^m*f.

  procedure Evaluate_Derivatives
              ( a : in out DoblDobl_Complex_Matrices.Matrix;
                r : in natural32; c : in out natural32; nq,nv,k : in natural32;
                f : in Poly_Sys; z : in DoblDobl_Complex_Vectors.Vector );
  procedure Evaluate_Derivatives
              ( file : in file_type;
                a : in out DoblDobl_Complex_Matrices.Matrix;
                r : in natural32; c : in out natural32; nq,nv,k : in natural32;
                f : in Poly_Sys; z : in DoblDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Evaluates all derivatives up to order k of the system f at z,
  --   in the rows of the matrix a, starting at given row r and column c.
  --   With a file as extra argument, diagnostics are written.

  procedure Compute_Derivatives
              ( a : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                r : in natural32; c : in out natural32; nq,nv,k : in natural32;
                f : in Poly_Sys );
  procedure Compute_Derivatives
              ( file : in file_type;
                a : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                r : in natural32; c : in out natural32; nq,nv,k : in natural32;
                f : in Poly_Sys );

  -- DESCRIPTION :
  --   Computes all derivatives up to order k of the system f,
  --   in the rows of the matrix a, starting at given row r and column c.
  --   With a file as extra argument, diagnostics are written.

  procedure Evaluate_All_Derivatives
              ( a : in out DoblDobl_Complex_Matrices.Matrix;
                r,c : in natural32; nq,nv,k : in natural32;
                f : in Poly_Sys; z : in DoblDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Evaluates all derivatives of f multiplied with all monomials of
  --   degree k-1, storing in the matrix a starting at row r and column c.

  procedure Compute_All_Derivatives
              ( a : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                r,c : in natural32; nq,nv,k : in natural32; f : in Poly_Sys );

  -- DESCRIPTION :
  --   Computes all derivatives of f multiplied with all monomials of
  --   degree k-1, storing in the matrix a starting at row r and column c.

  procedure Evaluate_Highest_Order
              ( a : in out DoblDobl_Complex_Matrices.Matrix;
                r,c : in natural32; nq,nv,k : in natural32;
                f : in Poly_Sys; z : in DoblDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Evaluates the k-th order derivative of f at all monomial multiples of
  --   degree 1 to k-1, storing in a, starting at row r and column c.

  procedure Compute_Highest_Order
              ( a : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                r,c : in natural32; nq,nv,k : in natural32; f : in Poly_Sys );

  -- DESCRIPTION :
  --   Computes the k-th order derivative of f at all monomial multiples of
  --   degree 1 to k-1, storing in a, starting at row r and column c.

  procedure Evaluate_Monomial_Multiples
              ( a : in out DoblDobl_Complex_Matrices.Matrix;
                r,c,nq,nv,k,nc1 : in natural32;
                f : in Poly_Sys; z : in DoblDobl_Complex_Vectors.Vector );
  procedure Evaluate_Monomial_Multiples
              ( file : in file_type;
                a : in out DoblDobl_Complex_Matrices.Matrix;
                r,c,nq,nv,k,nc1 : in natural32;
                f : in Poly_Sys; z : in DoblDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Evaluates monomial multiples of the degree k with the system f
  --   in the rows of the matrix a, starting at given row r and column c.
  --   Write intermediate output to file, if file given as argument.

  procedure Compute_Monomial_Multiples
              ( a : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                r,c,nq,nv,k,nc1 : in natural32; f : in Poly_Sys );
  procedure Compute_Monomial_Multiples
              ( file : in file_type;
                a : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                r,c,nq,nv,k,nc1 : in natural32; f : in Poly_Sys );

  -- DESCRIPTION :
  --   Evaluates monomial multiples of the degree k with the system f
  --   in the rows of the matrix a, starting at given row r and column c.
  --   Write intermediate output to file, if file given as argument.

end DoblDobl_Nullity_Polynomials;
