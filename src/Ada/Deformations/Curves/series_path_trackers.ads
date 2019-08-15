with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Homotopy_Continuation_Parameters;

package Series_Path_Trackers is

-- DESCRIPTION :
--   A series path tracker runs Newton's method on power series to
--   compute Pade approximants to predict solutions in the path tracker.
--   The procedures in this package give access to such trackers
--   in double, double double, and quad double precision.

  procedure Standard_Write 
              ( file : in file_type; nq,nv : in natural32;
                idxpar : in integer32;
                sols : in Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters );
  procedure DoblDobl_Write 
              ( file : in file_type; nq,nv : in natural32;
                idxpar : in integer32;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters );
  procedure QuadDobl_Write 
              ( file : in file_type; nq,nv : in natural32;
                idxpar : in integer32;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters );

  -- DESCRIPTION :
  --   Writes target system, start system (if idxpar > 0), start solutions,
  --   and the homotopy continuation parameters in pars to file.

  -- REQUIRED :
  --   Standard_Homotopy, DoblDobl_Homotopy, or QuadDobl_Homotopy 
  --   is well defined.

  -- ON ENTRY :
  --   file     file, opened for output;
  --   nq       number of equations;
  --   nv       number of variables;
  --   idxpar   index of the continuation parameter, 0 if artificial;
  --   sols     start solutions;
  --   pars     values of the homotopy continuation parameters.

  procedure Standard_Run
              ( nq,nvr,idxpar : in integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                sols : in out Standard_Complex_Solutions.Solution_List;
                vrb : in integer32 := 0 );
  procedure DoblDobl_Run
              ( nq,nvr,idxpar : in integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                vrb : in integer32 := 0 );
  procedure QuadDobl_Run
              ( nq,nvr,idxpar : in integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   With a homotopy defined, runs the path tracker starting at the
  --   given solutions, in double, double double or quad double precision.

  -- ON ENTRY :
  --   nq       number of equations in the homotopy;
  --   nvr      number of variables in the homotopy;
  --   idxpar   index of the parameter in a natural parameter homotopy,
  --            in 1..nvr, or else 0 for an artificial parameter homotopy;
  --   pars     values of the homotopy continuation parameters;
  --   sols     start solutions;
  --   verbose  the verbose level.

  -- ON RETURN :
  --   sols     solutions at the end of the path.

  function Prompt_for_Artificial return boolean;

  -- DESCRIPTION :
  --   Asks the user whether the homotopy is an artificial parameter
  --   homotopy and return true if so, otherwise false is returned.

  function Prompt_for_Homogenization ( dim : natural32 ) return natural32;

  -- DESCRIPTION :
  --   Asks the user whether homogeneous coordinates need to be used.
  --   On return is a natural number.  There are three possibilities:
  --   1) 0 : for affine coordinates
  --   2) 1 : in 1-homogeneous coordinates, in ordinary projective space,
  --   3) 2 or higher : in multi-projective coordinates,
  --   in a multi-projective space, defined by a partition of the variables.
  --   The total number of variables is given in the value of dim.

  procedure Standard_Main ( verbose : in integer32 := 0 );
  procedure DoblDobl_Main ( verbose : in integer32 := 0 );
  procedure QuadDobl_Main ( verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Prompts the user for a homotopy, runs the series path trackers
  --   in standard double, double double, or quad double precision.

end Series_Path_Trackers;
