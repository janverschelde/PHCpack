with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Black_Box_Root_Refiners is

-- DESCRIPTION :
--   Wraps the root refiners with basic settins for the tolerances.

  procedure Refine_Roots
               ( file : in file_type;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in out Standard_Complex_Solutions.Solution_List );
  procedure Refine_Roots
               ( file : in file_type;
                 p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure Refine_Roots
               ( file : in file_type;
                 p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Writes the result of the root refining on the system p,
  --   starting at the roots in sols to file.
  --   Computations occur in double, double double, or quad double precision.

  procedure Refine_Roots
               ( file : in file_type;
                 p : in Standard_Complex_Laur_Systems.Laur_Sys;
                 sols : in out Standard_Complex_Solutions.Solution_List );
  procedure Refine_Roots
               ( file : in file_type;
                 p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 sols : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure Refine_Roots
               ( file : in file_type;
                 p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 sols : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Applies the root refinement for Laurent systems,
  --   which is more basic as no deflation is available.

  procedure Standard_Main ( infilename,outfilename : in string;
                            verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Defines phc -b -v, to refine solutions in double precision,
  --   without prompting for root refining parameters.
  --   The first two arguments are the names of the input and output file.
  --   The last argument is the verbose level.

  procedure DoblDobl_Main ( infilename,outfilename : in string;
                            verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Defines phc -b -v, to refine solutions in double double precision,
  --   without prompting for root refining parameters.
  --   The first two arguments are the names of the input and output file.
  --   The last argument is the verbose level.

  procedure QuadDobl_Main ( infilename,outfilename : in string;
                            verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Defines phc -b -v, to refine solutions in quad double precision,
  --   without prompting for root refining parameters.
  --   The first two arguments are the names of the input and output file.
  --   The last argument is the verbose level.

end Black_Box_Root_Refiners;
