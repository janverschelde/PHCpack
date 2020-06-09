with text_io;                           use text_io;
with String_Splitters;                  use String_Splitters;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Multprec_Complex_Poly_Systems;     use Multprec_Complex_Poly_Systems;
with Multprec_Complex_Laur_Systems;     use Multprec_Complex_Laur_Systems;
with Multprec_Complex_Solutions;        use Multprec_Complex_Solutions;

package Multprec_System_and_Solutions_io is

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
  --   Writes the system p and the solution list to file,
  --   monomials in the system are on separate lines.

end Multprec_System_and_Solutions_io;
