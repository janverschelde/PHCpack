with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package Standard_Linear_Poly_Solvers is

-- DESCRIPTION :
--   This package parses a polynomial system into a linear system
--   and solves it.  This is useful to the blackbox solver in PHCpack.

  function Is_Linear ( p : Poly_Sys ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the degree of every polynomial in p is one.

  procedure Coefficients ( p : in Poly_Sys; 
                           A : out Matrix; b : out Vector );

  -- DESCRIPTION :
  --   Parses the linear polynomial system into the format A*x = b.

  -- REQUIRED : Is_Linear(p).

  -- ON ENTRY :
  --   p         a linear polynomial system.

  -- ON RETURN :
  --   A         coefficient matrix, A'range(1) = p'range
  --             and A'range(2) = 1..#variables in p;
  --   b         vector with constant coefficients of the polynomials,
  --             b'range = p'range.

  function Square_Solve ( A : Matrix; b : Vector ) return Solution;

  -- DESCRIPTION :
  --   Solves a square linear system and returns its solution.

  procedure Solve ( p : in Poly_Sys; s : out Solution; 
                    fail : out boolean );

  -- DESCRIPTION :
  --   Attempts to solve a system which is suppose to be linear.

  -- ON ENTRY :
  --   p         a polynomial system.

  -- ON RETURN :
  --   s         the solution to p if p is linear and square;
  --   fail      true if p is nonlinear or nonsquare.

end Standard_Linear_Poly_Solvers;
