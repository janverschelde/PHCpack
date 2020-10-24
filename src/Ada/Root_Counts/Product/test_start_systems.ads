with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;

package Test_Start_Systems is

-- DESCRIPTION :
--   Test on the construction of total-degree based start systems.

  procedure Enumerate_and_Solve
              ( d : in Standard_Natural_Vectors.Vector;
                c : in Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Lexicographically enumerates all solutions of the total degree
  --   start system with degrees in d and coefficients in c.

  procedure Create_and_Solve_Start_System ( n : in natural32 );

  -- DESCRIPTION :
  --   Prompts for the n degrees of the polynomials in a system.
  --   Applies the enumerator to the solutions of a start system
  --   with the given degrees.

  procedure Create_and_Solve_Start_System ( p : in Poly_Sys );

  -- DESCRIPTION :
  --   Makes a start system for p and checks the residuals
  --   of the computed solutions of the start system.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Start_Systems;
