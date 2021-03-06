with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with TripDobl_Complex_VecVecs;
with TripDobl_Speelpenning_Convolutions;

package Test_TripDobl_Fabry is

-- DESCRIPTION :
--   Given a polynomial system, adds a parameter t to every coefficient
--   and runs the Newton's method on the power series.
--   The smallest ratio of the coefficients of the series will give
--   the convergence radius and the location for the nearest singularity.
--   The tests are performed in triple double precision.

  procedure TripDobl_Newton_Steps
              ( s : in TripDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in TripDobl_Complex_VecVecs.VecVec;
                dim,deg : in integer32; maxit : in natural32 := 100 );

  -- DESCRIPTION :
  --   Applies several Newton steps on the system of convolution circuits s,
  --   departing from the series coefficients in scf.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system with solutions
  --   and defines the setup.

end Test_TripDobl_Fabry;
