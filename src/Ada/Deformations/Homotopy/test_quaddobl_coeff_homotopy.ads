with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Poly_Systems;

package Test_QuadDobl_Coeff_Homotopy is

-- DESCRIPTION :
--   Development of the evaluation of (1-t)*f + t*g
--   where monomials of f and g are shared
--   as in a coefficient-parameter homotopy, in quad double precision.

  procedure Write_Elements ( A,B : in QuadDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Writes the elements of the matrices A and B to screen,
  --   next to each other, for easy comparisons.

  procedure QuadDobl_Compared_Encapsulated_Eval ( n : in natural32 );

  -- DESCRIPTION :
  --   Evaluates (1-t)*p + t*q at some random t in [0,1]
  --   and at a random vector in quad double complex arithmetic.

  -- ON ENTRY :
  --   n        number of variables.

  procedure QuadDobl_Random_Systems
              ( n : in integer32;
                p,q : out QuadDobl_Complex_Poly_Systems.Poly_Sys );

  -- DESCRIPTION :
  --   Generates two n-dimensional systems p and q,
  --   with quad double complex coefficients.

  procedure QuadDobl_Random_Coefficient_Systems
              ( n : in integer32;
                p,q : out QuadDobl_Complex_Poly_Systems.Poly_Sys );

  -- DESCRIPTION :
  --   Generates two n-dimensional systems p and q,
  --   with the same supports.

  procedure QuadDobl_Compared_Encapsulation_Test;

  -- DESCRIPTION :
  --   Performs the evaluation test on randomly generated systems,
  --   comparing with the QuadDobl_Homotopy package.

  procedure QuadDobl_Homotopy_Performance ( n,m : natural32 );

  -- DESCRIPTION :
  --   Generates m random values for x and t
  --   for evaluation in the homotopy with timings taken.

  procedure QuadDobl_Coefficient_Homotopy_Performance ( n,m : natural32 );

  -- DESCRIPTION :
  --   Generates m random values for x and t
  --   for evaluation in the coefficient homotopy with timings taken.

  procedure QuadDobl_Performance_Test;

  -- DESCRIPTION :
  --   Performs the evaluation test on randomly generated systems,
  --   comparing with the QuadDobl_Homotopy package.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_QuadDobl_Coeff_Homotopy;
