with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_VecVecs;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with Homotopy_Continuation_Parameters;

package Multitasked_Path_Convolutions is

-- DESCRIPTION :
--   Multitasked path tracking with predictor-corrector-shift loops
--   on homotopy systems of convolution circuits.

  procedure Allocate ( v : in out Standard_Integer_VecVecs.VecVec;
                       n : in integer32 );

  -- DESCRIPTION :
  --   Allocates vectors of range 1..n in v.

  procedure Allocate ( v : in out Standard_Complex_VecVecs.VecVec;
                       n : in integer32 );
  procedure Allocate ( v : in out DoblDobl_Complex_VecVecs.VecVec;
                       n : in integer32 );
  procedure Allocate ( v : in out QuadDobl_Complex_VecVecs.VecVec;
                       n : in integer32 );

  -- DESCRIPTION :
  --   Allocates vectors of range 1..n in v,
  --   in double, double double, or quad double precision.

  procedure Standard_Multitasked_Tracker
              ( nbtasks : in integer32;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                hcrd,verbose : in boolean := true );
  procedure DoblDobl_Multitasked_Tracker
              ( nbtasks : in integer32;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                hcrd,verbose : in boolean := true );
  procedure QuadDobl_Multitasked_Tracker
              ( nbtasks : in integer32;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                hcrd,verbose : in boolean := true );

  -- DESCRIPTION :
  --   Applies multitasking to track all paths
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   nbtasks  the number of tasks;
  --   hom      system of homotopy convolution circuits;
  --   abh      radii as coefficients for mixed residuals;
  --   sols     start solutions;
  --   pars     values for the tolerances and parameters;
  --   hcrd     true if the homotopy is 1-homogeneous with one linear
  --            equation as its last equation which is then updated in
  --            the scaling of the solution coordinates,
  --            false if affine coordinates are used;
  --   verbose  indicates if extra output is requested.
  
  -- ON RETURN :
  --   sols     solutions at the end of the paths.

end Multitasked_Path_Convolutions;
