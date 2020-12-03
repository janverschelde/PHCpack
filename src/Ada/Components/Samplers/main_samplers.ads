with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Main_Samplers is

-- DESCRIPTION :
--   Wraps the samplers in local intrinsic coordinates,
--   in double, double double, and quad double precision.

  procedure Read_Witness_Set
              ( w : in string;
                p : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : out Standard_Complex_Solutions.Solution_List;
                dim : out natural32 );
  procedure Read_Witness_Set
              ( w : in string;
                p : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                dim : out natural32 );
  procedure Read_Witness_Set
              ( w : in string;
                p : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                dim : out natural32 );

  -- DESCRIPTION :
  --   Reads witness set k from file with name in the string w,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   w        name of file, if empty, the user is asked for a name.

  -- ON RETURN :
  --   p        polynomial equations defining the witness set;
  --   sols     points in the witness set;

  procedure Local_Intrinsic_Continuation
              ( file : in file_type; name : in string;
                n,k : in integer32; output : in boolean;
                f : in Standard_Complex_Poly_Systems.Poly_Sys;
                start : in Standard_Complex_Matrices.Matrix;
                esols : in Standard_Complex_Solutions.Solution_List );
  procedure Local_Intrinsic_Continuation
              ( file : in file_type; name : in string;
                n,k : in integer32; output : in boolean;
                f : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                start : in DoblDobl_Complex_Matrices.Matrix;
                esols : in DoblDobl_Complex_Solutions.Solution_List );
  procedure Local_Intrinsic_Continuation
              ( file : in file_type; name : in string;
                n,k : in integer32; output : in boolean;
                f : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                start : in QuadDobl_Complex_Matrices.Matrix;
                esols : in QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Deforms the given plane p into a random target plane,
  --   running intrinsic continuation on all solutions,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     file to write diagnostics on;
  --   name     name of the output file;
  --   n        ambient dimension, number of original variables;
  --   k        codimension of the solution set;
  --   output   if intermediate output wanted during tracking;
  --   f        original polynomial system;
  --   esols    extrinsic coordinates of the solutions.

  procedure Setup_Local_Coordinates
              ( file : in file_type; name : in string;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                d : in integer32 );
  procedure Setup_Local_Coordinates
              ( file : in file_type; name : in string;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                d : in integer32 );
  procedure Setup_Local_Coordinates
              ( file : in file_type; name : in string;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                d : in integer32 );

  -- DESCRIPTION :
  --   Prepares the setup for working with local intrinsic coordinates,
  --   in double, double double, or quad double precision.

  -- ON ENTRY:
  --   file     file to write diagnostics on;
  --   name     name of the output file;
  --   ep       embedded polynomial system with d slack variables;
  --   sols     generic points on the algebraic set;
  --   d        dimension of the witness set.

  procedure Sample_in_Standard_Precision ( witset,logfile : in string );
  procedure Sample_in_DoblDobl_Precision ( witset,logfile : in string );
  procedure Sample_in_QuadDobl_Precision ( witset,logfile : in string );

  -- DESCRIPTION :
  --   Prompts for a witness set and then computes
  --   a new set of samples using local intrinsic coordinates,
  --   for a new random set of slices in double, double double,
  --   or quad double precision.

  procedure Main ( witset,logfile : in string );

  -- DESCRIPTION :
  --   Given a witness set, computes new sample points on
  --   the algebraic set represented by the witness set.
  --   Defines phc -y.

  -- ON ENTRY :
  --   witset       witness set for an algebraic set;
  --   logfile      file name to write diagnostics on.

  -- ON RETURN :
  --   The output file contains new samples of an algebraic set,
  --   in witness set format.

end Main_Samplers;
