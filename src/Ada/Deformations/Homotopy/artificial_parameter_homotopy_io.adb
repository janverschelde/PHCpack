with text_io;                            use text_io;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with TripDobl_Complex_Poly_Systems_io;   use TripDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Standard_System_and_Solutions_io;
with DoblDobl_System_and_Solutions_io;
with TripDobl_System_and_Solutions_io;
with QuadDobl_System_and_Solutions_io;

package body Artificial_Parameter_Homotopy_io is

  procedure get ( target : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                  start : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                  sols : out Standard_Complex_Solutions.Solution_List ) is
  begin
    new_line;
    put_line("Reading the target system ..."); get(target);
    new_line;
    put_line("Reading the start system and its solutions ...");
    Standard_System_and_Solutions_io.get(start,sols);
  end get;

  procedure get ( target : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                  start : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                  sols : out DoblDobl_Complex_Solutions.Solution_List ) is
  begin
    new_line;
    put_line("Reading the target system ..."); get(target);
    new_line;
    put_line("Reading the start system and its solutions ...");
    DoblDobl_System_and_Solutions_io.get(start,sols);
  end get;

  procedure get ( target : out TripDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                  start : out TripDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                  sols : out TripDobl_Complex_Solutions.Solution_List ) is
  begin
    new_line;
    put_line("Reading the target system ..."); get(target);
    new_line;
    put_line("Reading the start system and its solutions ...");
    TripDobl_System_and_Solutions_io.get(start,sols);
  end get;

  procedure get ( target : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                  start : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                  sols : out QuadDobl_Complex_Solutions.Solution_List ) is
  begin
    new_line;
    put_line("Reading the target system ..."); get(target);
    new_line;
    put_line("Reading the start system and its solutions ...");
    QuadDobl_System_and_Solutions_io.get(start,sols);
  end get;

end Artificial_Parameter_Homotopy_io;
