with text_io;                            use text_io;
with String_Splitters;                   use String_Splitters;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Black_Box_Helpers is

-- DESCRIPTION :
--   Functions and procedures to help in the blackbox solvers.

  function Is_Constant_In 
              ( p : Standard_Complex_Polynomials.Poly ) return boolean;
  function Is_Constant_In 
              ( p : DoblDobl_Complex_Polynomials.Poly ) return boolean;
  function Is_Constant_In 
              ( p : QuadDobl_Complex_Polynomials.Poly ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the polynomial has a constant term.

  function Are_Constants_In
              ( p : Standard_Complex_Poly_Systems.Poly_Sys ) return boolean;
  function Are_Constants_In
              ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys ) return boolean;
  function Are_Constants_In
              ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys ) return boolean;

  -- DESCRIPTION :
  --   Returns true if all polynomials in p have a constant term.

  procedure Timing_Summary
              ( file : in file_type; roco,hoco,poco,total : in duration );

  -- DESCRIPTION :
  --   Writes a summary about execution times to the output file.

  procedure Append_Solutions_to_Input_File
              ( infilename : in string;
                sols : in Standard_Complex_Solutions.Solution_list;
                append_sols : in boolean );
  procedure Append_Solutions_to_Input_File
              ( infilename : in string;
                sols : in DoblDobl_Complex_Solutions.Solution_list;
                append_sols : in boolean );
  procedure Append_Solutions_to_Input_File
              ( infilename : in string;
                sols : in QuadDobl_Complex_Solutions.Solution_list;
                append_sols : in boolean );

  -- DESCRIPTION :
  --   If the solution list is not empty and append_sols is true,
  --   then the file with name in "infilename" is opened in append mode
  --   and the solutions are then appended to the input file.

  procedure Ask_Output_File
              ( outfile : out file_type; outfilename : in string;
                output_to_file : out boolean );
  procedure Ask_Output_File
              ( outfile : out file_type; outfilename : in string;
                output_to_file : out boolean;
                outnewname : out Link_to_String );

  -- DESCRIPTION :
  --   In case the output file is empty, the user is asked whether
  --   the output should be written to file.  In case the solutions
  --   should be written to file, the file with the given name is
  --   created, eventually after asking for a nonempty string.
  --   The new name of the string is optionally returned in outnewname.
  --   On return output_to_file is true if a file has been created.

end Black_Box_Helpers;
