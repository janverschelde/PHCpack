with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Vectors;
with Standard_Integer_Matrices;
with Standard_Complex_Matrices;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;

package Standard_Simplex_Systems is

-- DESCRIPTION :
--   We call a polynomial system a simplex system if by row reduction on 
--   the coefficient matrix, eventually after monomial divisions, the system
--   can be reduced to a binomial system.  A simplex system is represented
--   as C*x^A - b = 0, where A is an integer matrix which contains in its
--   columns the exponent vectors of the unknowns in x, where C is a general
--   matrix of complex coefficients and b is the constant righthand side.

-- FORMAT of a SIMPLEX SYSTEM : p(x) = 0 => C*x^A = b

  procedure Parse ( p : in Poly_Sys; nv : in integer32;
                    A : out Standard_Integer_Matrices.Matrix;
                    C : out Standard_Complex_Matrices.Matrix;
                    b : out Standard_Complex_Vectors.Vector;
                    fail : out boolean );

  procedure Parse ( p : in Laur_Sys; nv : in integer32;
                    A : out Standard_Integer_Matrices.Matrix;
                    C : out Standard_Complex_Matrices.Matrix;
                    b : out Standard_Complex_Vectors.Vector;
                    fail : out boolean );

  -- DESCRIPTION :
  --   Parses the equations of p into the format C*x^A = b.

  -- REQUIRED :
  --   p'range = A'range(2) = C'range(1) is 1..nq, nq = #equations;
  --   A'range(1) = C'range(2) is 1..nv, nv = #variables.

  -- ON ENTRY :
  --   p        a polynomial system in nv unknowns;
  --   nv       number of variables in p.

  -- ON RETURN :
  --   A        exponent vectors of x in its columns, if not fail;
  --   C        coefficient matrix of the system, if not fail;
  --   b        constant righthand side vector, if not fail;
  --   fail     true if the system has more than the allowed #monomials
  --            for it to be a simplex system.

  function Create ( A : Standard_Integer_Matrices.Matrix;
                    C : Standard_Complex_Matrices.Matrix;
                    b : Standard_Complex_Vectors.Vector ) return Poly_Sys;

  function Create ( A : Standard_Integer_Matrices.Matrix;
                    C : Standard_Complex_Matrices.Matrix;
                    b : Standard_Complex_Vectors.Vector ) return Laur_Sys;

  -- DESCRIPTION :
  --   Returns the system p(x) = C*x^A - b = 0, eventually after
  --   monomial multiplications to avoid negative exponents.

-- EVALUATION of a SIMPLEX SYSTEM :

  function Eval ( A : Standard_Integer_Matrices.Matrix;
                  C : Standard_Complex_Matrices.Matrix;
                  b,x : Standard_Complex_Vectors.Vector )
                return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the value of C*x^A - b.

end Standard_Simplex_Systems;
