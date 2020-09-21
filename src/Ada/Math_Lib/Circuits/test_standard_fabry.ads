with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_VecVecs;
with Standard_Complex_VecVecs;
with Standard_Speelpenning_Convolutions;
with Standard_Coefficient_Convolutions;

package Test_Standard_Fabry is

-- DESCRIPTION :
--   Given a polynomial system, adds a parameter t to every coefficient
--   and runs the Newton's method on the power series.
--   The smallest ratio of the coefficients of the series will give
--   the convergence radius and the location for the nearest singularity.
--   The tests are performed in standard double precision.

  procedure Standard_Newton_Steps
              ( s : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                dim,deg : in integer32; maxit : in natural32 := 100 );

  -- DESCRIPTION :
  --   Applies several Newton steps on the system of coefficient convolution
  --   circuits s, departing from the series coefficients in scf.

  procedure Standard_Newton_Steps
              ( s : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                dim,deg : in integer32; maxit : in natural32 := 100 );

  -- DESCRIPTION :
  --   Applies several Newton steps on the system of convolution circuits s,
  --   departing from the series coefficients in scf.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for a polynomial system with solutions, defines
  --   the setup and asks to run on coefficient convolutions or not.

end Test_Standard_Fabry;
