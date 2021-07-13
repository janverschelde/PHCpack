with text_io;                           use text_io;
with String_Splitters;                  use String_Splitters;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;

package Standard_System_and_Solutions_io is

-- DESCRIPTION :
--   A polynomial system is often accompanied by a list of solutions.
--   The routines in this package contain utilities to read and write
--   the system with its solution from and to file.

  procedure get ( n,m : out natural32; p : out Link_to_Array_of_Strings;
                  sols : out Solution_List;
                  banner : in string := "SOLUTIONS" );
  procedure get ( file : in file_type;
                  n,m : out natural32; p : out Link_to_Array_of_Strings;
                  sols : out Solution_List;
                  banner : in string := "SOLUTIONS" );
  procedure get ( p : out Link_to_Poly_Sys; sols : out Solution_List;
                  banner : in string := "SOLUTIONS" );
  procedure get ( file : in file_type;
                  p : out Link_to_Poly_Sys; sols : out Solution_List;
                  banner : in string := "SOLUTIONS" );
  procedure get ( p : out Link_to_Laur_Sys; sols : out Solution_List;
                  banner : in string := "SOLUTIONS" );
  procedure get ( file : in file_type;
                  p : out Link_to_Laur_Sys; sols : out Solution_List;
                  banner : in string := "SOLUTIONS" );

  -- DESCRIPTION :
  --   Reads a polynomial system and its solutions.
  --   By default, reading of solutions starts after the line
  --   that contains the banner string.
  --   Alternatives to the default banner are "START SOLUTIONS",
  --   to read the start solutions, or "THE SOLUTIONS",
  --   for the solutions of the system in case the file contains
  --   both a start system with start solutions and the target system
  --   and its solutions.

  procedure put ( file : in file_type;
                  p : in Poly_Sys; sols : in Solution_List;
                  banner : in string := "THE SOLUTIONS :" );
  procedure put ( file : in file_type;
                  p : in Laur_Sys; sols : in Solution_List;
                  banner : in string := "THE SOLUTIONS :" );

  -- DESCRIPTION :
  --   Writes the system p to file and if the solution list sols is not
  --   empty, then the solutions are written to the same file, preceded
  --   by the appropriate banner, recognizable by the get.
  --   The default banner assumes there is only one list of solutions.

  procedure put_line ( file : in file_type;
                       p : in Poly_Sys; sols : in Solution_List;
                       banner : in string := "THE SOLUTIONS :" );
  procedure put_line ( file : in file_type;
                       p : in Laur_Sys; sols : in Solution_List;
                       banner : in string := "THE SOLUTIONS :" );

  -- DESCRIPTION :
  --   Same as the put, but now the monomials in the system are
  --   each written to a separate line.

  procedure Scan_for_Start_System 
              ( infile : file_type; name : in Link_to_String;
                q : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                qsols : out Standard_Complex_Solutions.Solution_List;
                found : out boolean; verbose : in boolean := true;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Scans the infile for a start system and start solutions.

  -- ON ENTRY :
  --   infile   must be opened for input;
  --   name     name of the infile for error messages;
  --   verbose  if true, then error messages are displayed
  --            when a start system or its solutions are not found;
  --   vrblvl   is the verbose level, if positive then an opening
  --            message is written to screen.

  -- ON RETURN :
  --   q        a start system, if found;
  --   qsols    start solutions, if found;
  --   found    if true, then both q and qsols are found,
  --            otherwise, either q and/or qsols were not present.

  procedure Write_Scanned_Start_System
              ( name : in Link_to_String;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List );
  procedure Write_Scanned_Start_System
              ( file : in file_type;
                name : in Link_to_String;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Writes the scanned start system to standard output or to file.

  -- ON ENTRY :
  --   file     must be opened for output, or if omitted,
  --            the output will be sent to standard output;
  --   name     name of the file where the start system was scanned from;
  --   p        start system;
  --   sols     start solutions.

  procedure Main ( infilename,outfilename : in string;
                   vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Defines the phc -z -p option to extract the start system
  --   and its start solutions from the output of the blackbox solver.
  --   The names of the input and output files are given
  --   respectively in infilename and outfilename.

end Standard_System_and_Solutions_io;
