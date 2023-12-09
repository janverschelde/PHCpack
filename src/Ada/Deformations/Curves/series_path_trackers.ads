with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with TripDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with PentDobl_Complex_Solutions;
with OctoDobl_Complex_Solutions;
with DecaDobl_Complex_Solutions;
with HexaDobl_Complex_Solutions;
with Partitions_of_Sets_of_Unknowns;     use Partitions_of_Sets_of_Unknowns;
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
              ( monitor,verbose : in boolean;
                nq,nvr,idxpar : in integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                sols : in out Standard_Complex_Solutions.Solution_List;
                vrb : in integer32 := 0 );
  procedure DoblDobl_Run
              ( monitor,verbose : in boolean;
                nq,nvr,idxpar : in integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                vrb : in integer32 := 0 );
  procedure QuadDobl_Run
              ( monitor,verbose : in boolean;
                nq,nvr,idxpar : in integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   With a homotopy defined, runs the path tracker starting at the
  --   given solutions, in double, double double or quad double precision.
  --   No output is written to file.

  -- ON ENTRY :
  --   monitor  if a message is written at the start of each path;
  --   verbose  the verbose flag for extra output;
  --   nq       number of equations in the homotopy;
  --   nvr      number of variables in the homotopy;
  --   idxpar   index of the parameter in a natural parameter homotopy,
  --            in 1..nvr, or else 0 for an artificial parameter homotopy;
  --   pars     values of the homotopy continuation parameters;
  --   mhom     0 for affine, 1 for 1-homogenization, m for m-homogenization;
  --   idz      index representation of the partition z, for mhom > 1;
  --   sols     start solutions;
  --   vrb      the verbose level.

  -- ON RETURN :
  --   sols     solutions at the end of the path.

  procedure Standard_Run
              ( file : in file_type; monitor,verbose : in boolean;
                nq,nvr,idxpar : in integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                sols : in out Standard_Complex_Solutions.Solution_List;
                vrb : in integer32 := 0 );
  procedure DoblDobl_Run
              ( file : in file_type; monitor,verbose : in boolean;
                nq,nvr,idxpar : in integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                vrb : in integer32 := 0 );
  procedure QuadDobl_Run
              ( file : in file_type; monitor,verbose : in boolean;
                nq,nvr,idxpar : in integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   With a homotopy defined, runs the path tracker starting at the
  --   given solutions, in double, double double or quad double precision.

  -- ON ENTRY :
  --   file     file opened for output to write results to;
  --   monitor  if a message is written at the start of each path;
  --   verbose  the verbose flag for extra output;
  --   nq       number of equations in the homotopy;
  --   nvr      number of variables in the homotopy;
  --   idxpar   index of the parameter in a natural parameter homotopy,
  --            in 1..nvr, or else 0 for an artificial parameter homotopy;
  --   pars     values of the homotopy continuation parameters;
  --   mhom     0 for affine, 1 for 1-homogenization, m for m-homogenization;
  --   idz      index representation of the partition z, for mhom > 1;
  --   sols     start solutions;
  --   vrb      the verbose level.

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

  function Prompt_for_Partition
             ( nvr,mhom : in natural32 ) 
             return Standard_Natural_Vectors.Vector;

  -- DESCRIPTION :
  --   Prompts the user to define a partition of nvr number of variables.
  --   Returns a vector of range 1..nvr, defined as follows:
  --   at position k holds the index of set to which variable k belongs.
  --   This vector is the index representation of the partition.

  procedure Define_Partition
              ( n : in natural32; m : in out natural32;
                idx : out Standard_Natural_Vectors.Link_to_Vector;
                z : out Link_to_Partition );

  -- DESCRIPTION :
  --   Prompts the user for a partition of the set of n into m sets.
  --   The user is allowed to reset m, although m must be larger than one.
  --   On return in idx is the index representation of the partition z.

  procedure Write_Partition
              ( file : in file_type; n,m : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Writes the partition z of the set of n variables 
  --   for an m-homogenization to the output file.

  procedure Add_Multihomogeneous_Symbols
              ( m : in natural32; prefix : in string := "Z" );

  -- DESCRIPTION :
  --   Adds m symbols to represents the extra m symbols.
  --   If m = 1, then Z0 is added, otherwise, the m symbols
  --   start with 'Z' and have indices added, starting at 1, up to m.

  procedure Standard_Define_Homotopy
              ( nbq,nvr : out integer32;
                gamma : in Standard_Complex_Numbers.Complex_Number;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out Standard_Complex_Solutions.Solution_List );
  procedure DoblDobl_Define_Homotopy
              ( nbq,nvr : out integer32;
                gamma : in Standard_Complex_Numbers.Complex_Number;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure TripDobl_Define_Homotopy
              ( nbq,nvr : out integer32;
                gamma : in Standard_Complex_Numbers.Complex_Number;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out TripDobl_Complex_Solutions.Solution_List );
  procedure QuadDobl_Define_Homotopy
              ( nbq,nvr : out integer32;
                gamma : in Standard_Complex_Numbers.Complex_Number;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out QuadDobl_Complex_Solutions.Solution_List );
  procedure PentDobl_Define_Homotopy
              ( nbq,nvr : out integer32;
                gamma : in Standard_Complex_Numbers.Complex_Number;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out PentDobl_Complex_Solutions.Solution_List );
  procedure OctoDobl_Define_Homotopy
              ( nbq,nvr : out integer32;
                gamma : in Standard_Complex_Numbers.Complex_Number;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out OctoDobl_Complex_Solutions.Solution_List );
  procedure DecaDobl_Define_Homotopy
              ( nbq,nvr : out integer32;
                gamma : in Standard_Complex_Numbers.Complex_Number;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out DecaDobl_Complex_Solutions.Solution_List );
  procedure HexaDobl_Define_Homotopy
              ( nbq,nvr : out integer32;
                gamma : in Standard_Complex_Numbers.Complex_Number;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out HexaDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Prompts the user for target and start system and defines
  --   the artificial-parameter homotopy, eventually after applying
  --   projective transformations.

  -- ON ENTRY :
  --   gamma    gamma constant in the artificial-parameter homotopy.

  -- ON RETURN :
  --   nbq      number of equations;
  --   nvr      number of variables;
  --   mhom     0 for affine, 1 for 1-homogenization, m for m-homogenization;
  --   z        partition of the sets of unknowns, for mhom > 1;
  --   idz      index representation of the partition z, for mhom > 1;
  --   sols     start solutions.

  procedure Standard_Define_Homotopy
              ( nbq,nvr : out integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out Standard_Complex_Solutions.Solution_List );
  procedure DoblDobl_Define_Homotopy
              ( nbq,nvr : out integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure TripDobl_Define_Homotopy
              ( nbq,nvr : out integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out TripDobl_Complex_Solutions.Solution_List );
  procedure QuadDobl_Define_Homotopy
              ( nbq,nvr : out integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out QuadDobl_Complex_Solutions.Solution_List );
  procedure PentDobl_Define_Homotopy
              ( nbq,nvr : out integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out PentDobl_Complex_Solutions.Solution_List );
  procedure OctoDobl_Define_Homotopy
              ( nbq,nvr : out integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out OctoDobl_Complex_Solutions.Solution_List );
  procedure DecaDobl_Define_Homotopy
              ( nbq,nvr : out integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out DecaDobl_Complex_Solutions.Solution_List );
  procedure HexaDobl_Define_Homotopy
              ( nbq,nvr : out integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out HexaDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Prompts the user for target and start system and defines
  --   the artificial-parameter homotopy, eventually after applying
  --   projective transformations.

  -- ON ENTRY :
  --   pars     settings of the parameters (actually only for gamma).

  -- ON RETURN :
  --   nbq      number of equations;
  --   nvr      number of variables;
  --   mhom     0 for affine, 1 for 1-homogenization, m for m-homogenization;
  --   z        partition of the sets of unknowns, for mhom > 1;
  --   idz      index representation of the partition z, for mhom > 1;
  --   sols     start solutions.

  procedure Standard_Main ( vrb : in integer32 := 0 );
  procedure DoblDobl_Main ( vrb : in integer32 := 0 );
  procedure QuadDobl_Main ( vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Prompts the user for a homotopy, runs the series path trackers
  --   in standard double, double double, or quad double precision.

end Series_Path_Trackers;
