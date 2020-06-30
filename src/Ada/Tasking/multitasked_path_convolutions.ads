with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_VecVecs;
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
with Standard_Coefficient_Circuits;
with Standard_Coefficient_Convolutions;

package Multitasked_Path_Convolutions is

-- DESCRIPTION :
--   Multitasked path tracking with predictor-corrector-shift loops
--   on homotopy systems of convolution circuits.

  procedure Allocate ( v : in out Standard_Integer_VecVecs.VecVec;
                       n : in integer32 );

  -- DESCRIPTION :
  --   Allocates vectors of range 1..n in v.

  procedure Allocate ( v : in out Standard_Floating_VecVecs.VecVec;
                       n : in integer32 );
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
                hom : in Standard_Coefficient_Convolutions.Link_to_System;
                cfh,abh : in Standard_Coefficient_Circuits.Link_to_System;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Tracks all paths with coefficient convolution circuits.

  -- ON ENTRY :
  --   nbtasks  the number of tasks;
  --   hom      system of homotopy convolution circuits;
  --   cfh      circuits for the corrector;
  --   abh      radii as coefficients for mixed residuals;
  --   sols     start solutions;
  --   pars     values for the tolerances and parameters;
  --   mhom     0 if affine coordinates are used,
  --            1 for 1-homogeneous coordinates,
  --            m, for m > 1, for multi-homogenization;
  --   idz      the index representation of the partition of the variables,
  --            idz(k) returns a value between 1 and m,
  --            depending on which set the k-th variable belongs to;
  --   verbose  indicates if extra output is requested.
  
  -- ON RETURN :
  --   sols     solutions at the end of the paths.

  procedure Standard_Multitasked_Tracker
              ( nbtasks : in integer32;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                verbose : in boolean := true );
  procedure DoblDobl_Multitasked_Tracker
              ( nbtasks : in integer32;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                verbose : in boolean := true );
  procedure QuadDobl_Multitasked_Tracker
              ( nbtasks : in integer32;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Applies multitasking to track all paths
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   nbtasks  the number of tasks;
  --   hom      system of homotopy convolution circuits;
  --   abh      radii as coefficients for mixed residuals;
  --   sols     start solutions;
  --   pars     values for the tolerances and parameters;
  --   mhom     0 if affine coordinates are used,
  --            1 for 1-homogeneous coordinates,
  --            m, for m > 1, for multi-homogenization;
  --   idz      the index representation of the partition of the variables,
  --            idz(k) returns a value between 1 and m,
  --            depending on which set the k-th variable belongs to;
  --   verbose  indicates if extra output is requested.
  
  -- ON RETURN :
  --   sols     solutions at the end of the paths.

  procedure Track
              ( file : in file_type;
                hom : in Standard_Coefficient_Convolutions.Link_to_System;
                cfh,abh : in Standard_Coefficient_Circuits.Link_to_System;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbt,mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                arth : in boolean );

  -- DESCRIPTION :
  --   Tracks all paths starting at the solutions in sols,
  --   defined by the homotopy hom, with coefficient convolutions,
  --   in standard double precision.

  -- REQUIRED : the homotopy is square.

  -- ON ENTRY :
  --   file     file, opened for output;
  --   hom      system of homotopy convolution circuits;
  --   cfh      circuits for the corrector;
  --   abh      radii as coefficients for mixed residuals;
  --   sols     start solutions;
  --   pars     values for the tolerances and parameters;
  --   nbt      the number of tasks;
  --   mhom     0 if affine coordinates are used,
  --            1 for 1-homogeneous coordinates,
  --            m, for m > 1, for multi-homogenization;
  --   idz      the index representation of the partition of the variables,
  --            idz(k) returns a value between 1 and m,
  --            depending on which set the k-th variable belongs to;
  --   arth     true if the homotopy is an artificial-parameter one,
  --            false otherwise.

  -- ON RETURN :
  --   sols     solutions at the end of the path.

  procedure Track
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbt,mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                arth : in boolean );
  procedure Track
              ( file : in file_type;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbt,mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                arth : in boolean );
  procedure Track
              ( file : in file_type;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbt,mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                arth : in boolean );

  -- DESCRIPTION :
  --   Tracks all paths defined by the homotopy in hom,
  --   starting at solutions in sols,
  --   in double, double double, or quad double precision.

  -- REQUIRED : the homotopy is square.

  -- ON ENTRY :
  --   file     file, opened for output;
  --   hom      system of homotopy convolution circuits;
  --   abh      radii as coefficients for mixed residuals;
  --   sols     start solutions;
  --   pars     values for the tolerances and parameters;
  --   nbt      the number of tasks;
  --   mhom     0 if affine coordinates are used,
  --            1 for 1-homogeneous coordinates,
  --            m, for m > 1, for multi-homogenization;
  --   idz      the index representation of the partition of the variables,
  --            idz(k) returns a value between 1 and m,
  --            depending on which set the k-th variable belongs to;
  --   arth     true if the homotopy is an artificial-parameter one,
  --            false otherwise.

  -- ON RETURN :
  --   sols     solutions at the end of the path.

  procedure Standard_Main ( nt : in natural32; vrb : in integer32 := 0 );
  procedure DoblDobl_Main ( nt : in natural32; vrb : in integer32 := 0 );
  procedure QuadDobl_Main ( nt : in natural32; vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Prompts the user for a homotopy, setting of the parameters,
  --   tolerances, and the number of tasks, if nt = 0.
  --   If nt > 0, then nt will be the number of tasks.
  --   Then launches the trackers in standard double, double double, 
  --   or quad double precision.
  --   The verbose level is given in vrb.

end Multitasked_Path_Convolutions;
