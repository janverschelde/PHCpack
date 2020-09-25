with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with TripDobl_Complex_Matrices;
with TripDobl_Complex_Poly_Systems;

package Test_TripDobl_Coeff_Homotopy is

-- DESCRIPTION :
--   Development of the evaluation of (1-t)*f + t*g
--   where monomials of f and g are shared
--   as in a coefficient-parameter homotopy, in triple double precision.

  procedure Write_Elements ( A,B : in TripDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Writes the elements of the matrices A and B to screen,
  --   next to each other, for easy comparisons.

  procedure TripDobl_Compared_Encapsulated_Eval ( n : in natural32 );

  -- DESCRIPTION :
  --   Evaluates (1-t)*p + t*q at some random t in [0,1]
  --   and at a random vector in triple double complex arithmetic.

  -- ON ENTRY :
  --   n        number of variables.

  procedure TripDobl_Random_Systems
              ( n : in integer32;
                p,q : out TripDobl_Complex_Poly_Systems.Poly_Sys );

  -- DESCRIPTION :
  --   Generates two n-dimensional systems p and q,
  --   with triple double complex coefficients.

  procedure TripDobl_Random_Coefficient_Systems
              ( n : in integer32;
                p,q : out TripDobl_Complex_Poly_Systems.Poly_Sys );

  -- DESCRIPTION :
  --   Generates two n-dimensional systems p and q,
  --   with the same supports.

  procedure TripDobl_Compared_Encapsulation_Test;

  -- DESCRIPTION :
  --   Performs the evaluation test on randomly generated systems,
  --   comparing with the TripDobl_Homotopy package.

  procedure TripDobl_Homotopy_Performance ( n,m : natural32 );

  -- DESCRIPTION :
  --   Generates m random values for x and t
  --   for evaluation in the homotopy with timings taken.

  procedure TripDobl_Coefficient_Homotopy_Performance ( n,m : natural32 );

  -- DESCRIPTION :
  --   Generates m random values for x and t
  --   for evaluation in the coefficient homotopy with timings taken.

  procedure TripDobl_Performance_Test;

  -- DESCRIPTION :
  --   Performs the evaluation test on randomly generated systems,
  --   comparing with the TripDobl_Homotopy package.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_TripDobl_Coeff_Homotopy;
