with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with OctoDobl_Complex_Numbers;          use OctoDobl_Complex_Numbers;
with OctoDobl_Complex_Vectors;          use OctoDobl_Complex_Vectors;
with OctoDobl_Complex_Poly_Systems;     use OctoDobl_Complex_Poly_Systems;

package Test_OctoDobl_Poly_Homotopy is

-- DESCRIPTION :
--   Tests the polynomial homotopy in octo double precision.

  procedure Test_Homotopy ( p,q : in Poly_Sys );

  -- DESCRIPTION :
  --   Tests the memory management with frequent execution
  --   of the create and clear procedures,
  --   on the target p and start q.

  procedure Test_Evaluation
              ( p,q : in Poly_Sys; x : in Vector;
                t,gamma : in Complex_Number; k : in natural32 );

  -- DESCRIPTION :
  --   Evaluates the homotopy and compares with computed values.

  -- ON ENTRY :
  --   p        target system;
  --   q        start system;
  --   x        coordinates of some point;
  --   t        value of the continuation parameter;
  --   gamma    the gamma constant in the homotopy;
  --   k        power of the continuation parameter.

  procedure Random_Test_Homotopy_Eval ( p,q : in Poly_Sys );

  -- DESCRIPTION :
  --   Evaluates the homotopy at random points,
  --   for the homotopy with target p and start q.

  procedure Interactive_Test ( p,q : in Poly_Sys );

  -- DESCRIPTION :
  --   Prompts for complex numbers to evaluate the homotopy,
  --   for target p and start q.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for a target and start system,
  --   displays a menu and prompts for a test.

end Test_OctoDobl_Poly_Homotopy;
