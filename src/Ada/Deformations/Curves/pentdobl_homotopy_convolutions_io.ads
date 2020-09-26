with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Partitions_of_Sets_of_Unknowns;     use Partitions_of_Sets_of_Unknowns;
with PentDobl_Complex_Solutions;         use PentDobl_Complex_Solutions;
with PentDobl_Speelpenning_Convolutions; use PentDobl_Speelpenning_Convolutions;

package PentDobl_Homotopy_Convolutions_io is

-- DESCRIPTION :
--   Provides a general interactive procedure to prompt the user to
--   define a homotopy respresented by convolution circuits,
--   and corresponding start solutions, in penta double precision.

  function Make_Homotopy ( nq,idx,deg : integer32 ) return Link_to_System;

  -- DESCRIPTION :
  --   Returns a system of convolution circuits for the homotopy
  --   defined in penta double precision in PentDobl_Homotopy.

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

  procedure get ( deg : in integer32; artificial : in boolean;
                  gamma : in Complex_Number; h : out Link_to_System;
                  s : out Solution_List; idxpar : out integer32;
                  mhom : out natural32; z : out Link_to_Partition;
                  idz : out Standard_Natural_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Prompts the user for a homotopy for a given degree deg,
  --   and for the type of homogenization.
  --
  -- ON ENTRY :
  --   deg        degree of the power series coefficients in h;
  --   artificial is true if the homotopy is an artificial parameter
  --              homotopy, false otherwise.
  --
  -- ON RETURN :
  --   h          convolution circuits of degree deg;
  --   s          start solutions, in a natural parameter homotopy,
  --              the value for the continuation has been dropped
  --              from the solutions given on input;
  --   idxpar     the index of the homotopy continuation parameter,
  --              in a natural parameter homotopy, idxpar is in the
  --              range 1..h.dim, otherwise idxpar is h.neq + 1,
  --              assuming a square polynomial system;
  --   mhom       0 for affine, 1 for 1-homogenization,
  --              equals some m larger than 1, for m-homogenization;
  --   z          partition of the sets of unknowns, for mhom > 1;
  --   idz        index representation of the partition z, for mhom > 1.

end PentDobl_Homotopy_Convolutions_io;
