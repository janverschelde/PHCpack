with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with QuadDobl_Complex_Vectors;         use QuadDobl_Complex_Vectors;
with Standard_Integer64_Matrices;
with Multprec_Integer_Matrices;
with QuadDobl_Complex_Solutions;       use QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Poly_Systems;    use QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;    use QuadDobl_Complex_Laur_Systems;

package QuadDobl_Binomial_Systems is

-- DESCRIPTION :
--   A binomial system is a system with exactly two terms with nonzero
--   (standard complex) coefficient in every equation.
--   We represent a binomial system as x^A - c = 0, where the matrix A
--   contains in its columns the exponent vector of the unknowns in x.
--   This package offers functions to define and evaluate a binomial system.

-- FORMAT of a BINOMIAL SYSTEM :  p(x) = 0 => x^A = c

  procedure Parse ( p : in Poly_Sys; nq : in integer32;
                    A : out Standard_Integer64_Matrices.Matrix;
                    c : out Vector; fail : out boolean );

  procedure Parse ( p : in Laur_Sys; nq : in integer32;
                    A : out Standard_Integer64_Matrices.Matrix;
                    c : out Vector; fail : out boolean );

  -- DESCRIPTION :
  --   Parses the equations of p into the format x^A = c.

  -- REQUIRED :
  --   p'range = A'range(2) = c'range = 1..nq, nq = #equations;
  --   A'range(1) = 1..nv, nv = #equations.

  -- ON ENTRY :
  --   p        a polynomial system, not necessarily square;
  --   nq       number of equations in p.

  -- ON RETURN :
  --   A        the exponent vectors in the columns, if not fail;
  --   c        coefficient vector, if not fail;
  --   fail     true if not exactly two monomials in every equation,
  --            false otherwise.

  function Create ( A : Standard_Integer64_Matrices.Matrix;
                    c : Vector ) return Poly_Sys;
  function Create ( A : Standard_Integer64_Matrices.Matrix;
                    c : Vector ) return Laur_Sys;

  -- DESCRIPTION :
  --   Returns the system p(x) = x^A - c = 0,
  --   modulo monomial multiplications to avoid negative exponents.

  -- REQUIRED : A'range(2) = c'range.

-- EVALUATION of a BINOMIAL SYSTEM :

  function Eval ( A : Standard_Integer64_Matrices.Matrix;
                  x : Vector ) return Vector;
  function Eval ( A : Multprec_Integer_Matrices.Matrix;
                  x : Vector ) return Vector;

  -- DESCRIPTION : returns x^A.

  procedure Eval ( A : in Standard_Integer64_Matrices.Matrix;
                   x : in Vector; y : out Vector );

  -- DESCRIPTION :
  --   Returns in y = x^A.

  -- REQUIRED : y'range = A'range(2).

  function Eval ( A : Standard_Integer64_Matrices.Matrix;
                  c,x : Vector ) return Vector;

  -- DESCRIPTION : returns x^A - c.

  function Eval ( A : Standard_Integer64_Matrices.Matrix;
                  s : Solution_List ) return Solution_List;
  function Eval ( A : Multprec_Integer_Matrices.Matrix;
                  s : Solution_List ) return Solution_List;

  -- DESCRIPTION :
  --   The solutions on return are vectors Eval(A,x), for all x in s.
  --   For a multiprecision A, all components of x must have modulus 1.

  procedure Eval ( A : in Standard_Integer64_Matrices.Matrix;
                   s : in Solution_List; w : in out Vector );

  -- DESCRIPTION :
  --   Replaces the solutions x in sols by Eval(A,x), for all x in s,
  --   using w as a work vector for intermediate products.

end QuadDobl_Binomial_Systems;
