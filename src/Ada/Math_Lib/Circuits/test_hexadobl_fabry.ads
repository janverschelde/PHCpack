with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with HexaDobl_Complex_VecVecs;
with HexaDobl_Speelpenning_Convolutions;

package Test_HexaDobl_Fabry is

-- DESCRIPTION :
--   Given a polynomial system, adds a parameter t to every coefficient
--   and runs the Newton's method on the power series.
--   The smallest ratio of the coefficients of the series will give
--   the convergence radius and the location for the nearest singularity.
--   The tests are performed in hexa double precision.

  procedure HexaDobl_Newton_Steps
              ( s : in HexaDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in HexaDobl_Complex_VecVecs.VecVec;
                dim,deg : in integer32; maxit : in natural32 := 100 );

  -- DESCRIPTION :
  --   Applies several Newton steps on the system of convolution circuits s,
  --   departing from the series coefficients in scf.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for a polynomial system with solutions
  --   and defines the setup.

end Test_HexaDobl_Fabry;
