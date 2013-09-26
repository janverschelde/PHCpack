with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package Standard_Sparse_Solvers is

-- DESCRIPTION :
--   This package provides drivers to the fewnomial system solvers.

  function Is_Fewnomial_System ( p : Poly_Sys ) return boolean;
  function Is_Fewnomial_System ( p : Laur_Sys ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the given system is a fewnomial system,
  --   i.e.: it has fewer than n+1 distinct monomials,
  --   eventually after shifting monomials by division.

  procedure Solve ( p : in Poly_Sys; tol : in double_float;
                    sols : out Solution_List;
                    fail,zero_y : out boolean );
  procedure Solve ( p : in Laur_Sys; tol : in double_float;
                    sols : out Solution_List;
                    fail,zero_y : out boolean );

  procedure Solve ( p : in Poly_Sys; tol : in double_float;
                    rcond : out double_float; sols : out Solution_List;
                    fail,zero_y : out boolean );
  procedure Solve ( p : in Laur_Sys; tol : in double_float;
                    rcond : out double_float; sols : out Solution_List;
                    fail,zero_y : out boolean );

  procedure Solve ( p : in Poly_Sys; tol : in double_float;
                    rcond : out double_float; sols : out Solution_List;
                    fail,zero_y : out boolean; rsum : out double_float );
  procedure Solve ( p : in Laur_Sys; tol : in double_float;
                    rcond : out double_float; sols : out Solution_List;
                    fail,zero_y : out boolean; rsum : out double_float );

  -- DESCRIPTION :
  --   Returns all solutions of p if p is a fewnomial system,
  --   otherwise false is returned.

  -- ON ENTRY :
  --   p        a polynomial system;
  --   tol      tolerance to decide whether the modulus of a complex
  --            number is small enough to be considered as zero.

  -- ON RETURN : 
  --   sols     solutions to p if p is a fewnomial system;
  --   rcond    estimate for inverse of condition number of the
  --            condition number of p, only meaningful if not fail;
  --   fail     true if p is not a fewnomial system;
  --   zero_y   true if the intermediate linear system has a solution with
  --            at least one component equal to zero, which implies that
  --            the system can have no solutions with all components
  --            different from zero, otherwise zero_y is false on return;
  --   rsum     sum of all residuals of the solutions at p.

end Standard_Sparse_Solvers;
