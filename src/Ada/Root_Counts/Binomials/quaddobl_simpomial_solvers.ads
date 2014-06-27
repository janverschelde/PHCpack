with Quad_Double_Numbers;                use Quad_Double_Numbers;
with QuadDobl_Complex_Poly_Systems;      use QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;      use QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Solutions;         use QuadDobl_Complex_Solutions;

package QuadDobl_Simpomial_Solvers is

-- DESCRIPTION :
--   This package provides drivers to the simplex system solvers.

  function Is_Simplex_System ( p : Poly_Sys ) return boolean;
  function Is_Simplex_System ( p : Laur_Sys ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the given system is a simplex system,
  --   i.e.: it has fewer than n+1 distinct monomials,
  --   eventually after shifting monomials by division.

  procedure Solve ( p : in Poly_Sys; tol : in quad_double;
                    sols : out Solution_List;
                    fail,zero_y : out boolean;
                    multprec_hermite : in boolean := false );
  procedure Solve ( p : in Laur_Sys; tol : in quad_double;
                    sols : out Solution_List;
                    fail,zero_y : out boolean;
                    multprec_hermite : in boolean := false );

  procedure Solve ( p : in Poly_Sys; tol : in quad_double;
                    rcond : out quad_double; sols : out Solution_List;
                    fail,zero_y : out boolean;
                    multprec_hermite : in boolean := false );
  procedure Solve ( p : in Laur_Sys; tol : in quad_double;
                    rcond : out quad_double; sols : out Solution_List;
                    fail,zero_y : out boolean;
                    multprec_hermite : in boolean := false );

  procedure Solve ( p : in Poly_Sys; tol : in quad_double;
                    rcond : out quad_double; sols : out Solution_List;
                    fail,zero_y : out boolean; rsum : out quad_double;
                    multprec_hermite : in boolean := false );
  procedure Solve ( p : in Laur_Sys; tol : in quad_double;
                    rcond : out quad_double; sols : out Solution_List;
                    fail,zero_y : out boolean; rsum : out quad_double;
                    multprec_hermite : in boolean := false );

  -- DESCRIPTION :
  --   Returns all solutions of p if p is a simplex system,
  --   otherwise false is returned.

  -- ON ENTRY :
  --   p        a polynomial system;
  --   tol      tolerance to decide whether the modulus of a complex
  --            number is small enough to be considered as zero;
  --   multprec_hermite indicates whether multiprecision arithmetic
  --            will be used to compute the Hermite normal form.

  -- ON RETURN : 
  --   sols     solutions to p if p is a simplex system;
  --   rcond    estimate for inverse of condition number of the
  --            condition number of p, only meaningful if not fail;
  --   fail     true if p is not a simplex system;
  --   zero_y   true if the intermediate linear system has a solution with
  --            at least one component equal to zero, which implies that
  --            the system can have no solutions with all components
  --            different from zero, otherwise zero_y is false on return;
  --   rsum     sum of all residuals of the solutions at p.

end QuadDobl_Simpomial_Solvers;
