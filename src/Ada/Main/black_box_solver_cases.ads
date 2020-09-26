with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Black_Box_Solver_Cases is

-- DESCRIPTION :
--   Special cases of polynomial systems have only one single equation,
--   are linear, or are supported on a simplex.

  procedure Solve_for_Special_Cases
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                rc : out natural32;
                sols : out Standard_Complex_Solutions.Solution_List;
                fail : out boolean; verbose : in integer32 := 0 );
  procedure Solve_for_Special_Cases
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                rc : out natural32;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                fail : out boolean; verbose : in integer32 := 0 );
  procedure Solve_for_Special_Cases
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                rc : out natural32;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                fail : out boolean; verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Checks whether the system p is one of the three special cases:
  --   1. a polynomial in one variable;
  --   2. a linear system;
  --   3. a binomial system with a nonzero constant.
  --   In one of these three cases, the solver will return solutions
  --   and fail will be false.  Otherwise, fail is true.

  -- ON INPUT :
  --   p            a polynomial system;
  --   verbose      the verbose level.

  -- ON RETURN :
  --   rc           equals the number of solutions in sols or 0 if fail;
  --   fail         true if the system is not a special case,
  --                false otherwise;
  --   sols         solutions of p if not fail.

end Black_Box_Solver_Cases;
