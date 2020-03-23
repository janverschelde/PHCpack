with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Speelpenning_Convolutions; use Standard_Speelpenning_Convolutions;

package Standard_Homotopy_Convolutions_io is

-- DESCRIPTION :
--   Provides a general interactive procedure to prompt the user to
--   define a homotopy respresented by convolution circuits,
--   and corresponding start solutions, in hardward double precision.

  function Make_Homotopy ( nq,idx,deg : integer32 ) return Link_to_System;

  -- DESCRIPTION :
  --   Returns a system of convolution circuits for the homotopy
  --   defined in standard double precision in Standard_Homotopy.

  -- ON ENTRY :
  --   nq         number of equations in the homotopy;
  --   idx        index of the continuation parameter;
  --   deg        degree of the power series coefficients.

  procedure get ( deg : in integer32; h : out Link_to_System;
                  s : out Solution_List; idxpar : out integer32 );

  -- DESCRIPTION :
  --   Prompts the user for a homotopy for a given degree deg.
  --
  -- ON ENTRY :
  --   deg        degree of the power series coefficients in h.
  --
  -- ON RETURN :
  --   h          convolution circuits of degree deg;
  --   s          start solutions, in a natural parameter homotopy,
  --              the value for the continuation has been dropped
  --              from the solutions given on input;
  --   idxpar     the index of the homotopy continuation parameter,
  --              in a natural parameter homotopy, idxpar is in the
  --              range 1..h.dim, otherwise idxpar is h.neq + 1,
  --              assuming a square polynomial system.

end Standard_Homotopy_Convolutions_io;
