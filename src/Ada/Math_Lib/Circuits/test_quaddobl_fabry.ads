with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Speelpenning_Convolutions;

package Test_QuadDobl_Fabry is

-- DESCRIPTION :
--   Given a polynomial system, adds a parameter t to every coefficient
--   and runs the Newton's method on the power series.
--   The smallest ratio of the coefficients of the series will give
--   the convergence radius and the location for the nearest singularity.
--   The tests are performed in quad double precision.

  procedure QuadDobl_Newton_Steps
              ( s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in QuadDobl_Complex_VecVecs.VecVec;
                dim,deg : in integer32; maxit : in natural32 := 100 );

  -- DESCRIPTION :
  --   Applies several Newton steps on the system of convolution circuits s,
  --   departing from the series coefficients in scf.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system with solutions
  --   and defines the setup.

end Test_QuadDobl_Fabry;
