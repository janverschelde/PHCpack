with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with Homotopy_Continuation_Parameters;

package Track_Path_Convolutions is

-- DESCRIPTION :
--   This package provides the main interactive procedures to launch
--   the path trackers on homotopies defined by convolution circuits.

  procedure Track
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                arth : in boolean );
  procedure Track
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                arth : in boolean );
  procedure Track
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                arth : in boolean );

  -- DESCRIPTION :
  --   Tracks all paths defined by the homotopy in hom,
  --   starting at solutions in sols,
  --   in double, double double, or quad double precision.

  -- REQUIRED : the homotopy is square.

  -- ON ENTRY :
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
  --   arth     true if the homotopy is an artificial-parameter one,
  --            false otherwise.

  -- ON RETURN :
  --   sols     solutions at the end of the path.

  procedure Main
              ( hom : out Standard_Speelpenning_Convolutions.Link_to_System;
                abh : out Standard_Speelpenning_Convolutions.Link_to_System;
                arth : out boolean;
                pars : out Homotopy_Continuation_Parameters.Parameters;
                sols : out Standard_Complex_Solutions.Solution_List;
                mhom : out natural32;
                idz : out Standard_Natural_Vectors.Link_to_Vector );
  procedure Main
              ( hom : out DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : out DoblDobl_Speelpenning_Convolutions.Link_to_System;
                arth : out boolean;
                pars : out Homotopy_Continuation_Parameters.Parameters;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                mhom : out natural32;
                idz : out Standard_Natural_Vectors.Link_to_Vector );
  procedure Main
              ( hom : out QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : out QuadDobl_Speelpenning_Convolutions.Link_to_System;
                arth : out boolean;
                pars : out Homotopy_Continuation_Parameters.Parameters;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                mhom : out natural32;
                idz : out Standard_Natural_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Promps the user for a homotopy, tunes the parameter settings,
  --   in double, double double, or quad double precision.

  -- ON RETURN :
  --   hom      convolution circuits for a homotopy;
  --   abh      absolute value coefficient homotopy for mixed residuals;
  --   arth     flag to indicatei if homotopy is artificial parameter,
  --            or (if false), if hom is a natural parameter homotopy;
  --   pars     tuned values of the settings of parameters and tolerances;
  --   sols     start solutions for the homotopy;
  --   mhom     0 if affine, m if m-homogeneous, for m > 0;
  --   idz      index representation of the partition of the variables,
  --            only if mhom > 1, then idz(k) is the index of the set
  --            the k-th variable belongs to.

  procedure Standard_Main ( vrb : in integer32 := 0 );
  procedure DoblDobl_Main ( vrb : in integer32 := 0 );
  procedure QuadDobl_Main ( vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Prompts the user for a homotopy, then launches the trackers
  --   in standard double, double double, or quad double precision.
  --   The verbose level is given in vrb.

end Track_Path_Convolutions;
