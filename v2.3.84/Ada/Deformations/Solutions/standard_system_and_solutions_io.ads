with text_io;                           use text_io;
with String_Splitters;                  use String_Splitters;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;

package Standard_System_and_Solutions_io is

-- DESCRIPTION :
--   A polynomial system is often accompanied by a list of solutions.
--   The routines in this package contain utilities to read and write
--   the system with its solution from and to file.

  procedure get ( n,m : out natural32; p : out Link_to_Array_of_Strings;
                  sols : out Solution_List );
  procedure get ( file : in file_type;
                  n,m : out natural32; p : out Link_to_Array_of_Strings;
                  sols : out Solution_List );
  procedure get ( p : out Link_to_Poly_Sys; sols : out Solution_List );
  procedure get ( file : in file_type;
                  p : out Link_to_Poly_Sys; sols : out Solution_List );
  procedure get ( p : out Link_to_Laur_Sys; sols : out Solution_List );
  procedure get ( file : in file_type;
                  p : out Link_to_Laur_Sys; sols : out Solution_List );

  -- DESCRIPTION :
  --   Reads a polynomial system and its solutions.

  procedure put ( file : in file_type;
                  p : in Poly_Sys; sols : in Solution_List );
  procedure put ( file : in file_type;
                  p : in Laur_Sys; sols : in Solution_List );

  -- DESCRIPTION :
  --   Writes the system p to file and if the solution list sols is not
  --   empty, then the solutions are written to the same file, preceded
  --   by the appropriate banner, recognizable by the get.

  procedure put_line ( file : in file_type;
                       p : in Poly_Sys; sols : in Solution_List );
  procedure put_line ( file : in file_type;
                       p : in Laur_Sys; sols : in Solution_List );

  -- DESCRIPTION :
  --   Same as the put, but now the monomials in the system are
  --   each written to a separate line.

end Standard_System_and_Solutions_io;
