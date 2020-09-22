with text_io;                           use text_io;
with TripDobl_Complex_Poly_Systems;     use TripDobl_Complex_Poly_Systems;
with TripDobl_Complex_Laur_Systems;     use TripDobl_Complex_Laur_Systems;
with TripDobl_Complex_Solutions;        use TripDobl_Complex_Solutions;

package TripDobl_System_and_Solutions_io is

-- DESCRIPTION :
--   A polynomial system is often accompanied by a list of solutions.
--   The routines in this package contain utilities to write
--   the system with its solution to file.

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

end TripDobl_System_and_Solutions_io;
