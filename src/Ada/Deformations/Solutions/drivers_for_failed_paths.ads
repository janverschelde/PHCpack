with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Solutions;

package Drivers_for_Failed_Paths is

-- DESCRIPTION :
--   Provides tools to help processing the output files of trackers.

  procedure Prompt_for_File
              ( file : in out file_type; len : in natural32;
                tofile : out boolean );

  -- DESCRIPTION :
  --   Given the length of the solution list in len,
  --   prompts the user for a file name if len > 0
  --   and the created file object is returned in file.
  --   The file is opened for output if tofile is true on return.

  procedure Write_to_File
              ( file : in file_type; len : in natural32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List );
  procedure Write_to_File
              ( file : in file_type; len : in natural32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List );
  procedure Write_to_File
              ( file : in file_type; len : in natural32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List );
  procedure Write_to_File
              ( file : in file_type; len : in natural32;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List );
  procedure Write_to_File
              ( file : in file_type; len : in natural32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List );
  procedure Write_to_File
              ( file : in file_type; len : in natural32;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Writes the system and solutions to file.

  -- REQUIRED : len = Length_Of(sols) > 0.

  -- ON ENTRY :
  --   file     file opened for output;
  --   len      number of solutions in sols;
  --   p        a start system used in a homotopy.
  --   sols     solutions of p which correspond to failed paths.

  procedure Write_to_File
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List );
  procedure Write_to_File
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List );
  procedure Write_to_File
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List );
  procedure Write_to_File
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List );
  procedure Write_to_File
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List );
  procedure Write_to_File
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   If the solution list sols is not empty,
  --   then the user is prompted for a file name.
  --   Prompts the user for a file name to write the solutions on.
  --   If the length of sols is positive,
  --   then the system p and the solutions in sols are written to file.

  -- ON ENTRY :
  --   p        a start system used in a homotopy.
  --   sols     solutions of p which correspond to failed paths.

  procedure Standard_Scan_Failed_Paths
              ( infilename,outfilename : in string );
  procedure DoblDobl_Scan_Failed_Paths
              ( infilename,outfilename : in string );
  procedure QuadDobl_Scan_Failed_Paths
              ( infilename,outfilename : in string );

  -- DESCRIPTION :
  --   Scans the output file of a path tracker for failed paths,
  --   in double, double double, or quad double precision.
  --   If the names of the files are empty strings,
  --   then the user will be prompted for file names.

  -- ON ENTRY :
  --   infilename is the file name of the input file;
  --   outfilename is the file name of the output file.

end Drivers_for_Failed_Paths;
