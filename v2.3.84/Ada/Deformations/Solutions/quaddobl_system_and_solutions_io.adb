with Communications_with_User;           use Communications_with_User;
with File_Scanning;                      use File_Scanning;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Laur_Systems_io;   use QuadDobl_Complex_Laur_Systems_io;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;

package body QuadDobl_System_and_Solutions_io is

-- AUXILIARY ROUTINES :

  procedure Scan_for_Solutions
              ( file : in file_type; sols : out Solution_List ) is

    found : boolean;

  begin
    Scan_and_Skip(file,"SOLUTIONS",found);
    if found
     then get(file,sols);
    end if;
  end Scan_for_Solutions;

  procedure Write_Solutions
              ( file : in file_type; sols : in Solution_List ) is
  begin
    if not Is_Null(sols) then
      new_line(file);
      put_line(file,"THE SOLUTIONS : ");
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    end if;
  end Write_Solutions;

-- TARGET ROUTINES :

  procedure get ( p : out Link_to_Poly_Sys; sols : out Solution_List ) is

    file : file_type;

  begin
    new_line;
    put_line("Reading the name of the input file...");
    Read_Name_and_Open_File(file);
    get(file,p,sols);
    close(file);
  end get;

  procedure get ( file : in file_type;
                  p : out Link_to_Poly_Sys; sols : out Solution_List ) is
  begin
    get(file,p);
    Scan_for_Solutions(file,sols);
  end get;

  procedure get ( p : out Link_to_Laur_Sys; sols : out Solution_List ) is

    file : file_type;

  begin
    new_line;
    put_line("Reading the name of the input file...");
    Read_Name_and_Open_File(file);
    get(file,p,sols);
    close(file);
  end get;

  procedure get ( file : in file_type;
                  p : out Link_to_Laur_Sys; sols : out Solution_List ) is
  begin
    get(file,p);
    Scan_for_Solutions(file,sols);
  end get;

  procedure put ( file : in file_type;
                  p : in Poly_Sys; sols : in Solution_List ) is
  begin
    put(file,p'last,1); new_line(file);
    put(file,p);
    Write_Solutions(file,sols);
  end put;

  procedure put ( file : in file_type;
                  p : in Laur_Sys; sols : in Solution_List ) is
  begin
    put(file,p'last,1); new_line(file);
    put(file,p);
    Write_Solutions(file,sols);
  end put;

  procedure put_line ( file : in file_type;
                       p : in Poly_Sys; sols : in Solution_List ) is
  begin
    put_line(file,p);
    Write_Solutions(file,sols);
  end put_line;

  procedure put_line ( file : in file_type;
                       p : in Laur_Sys; sols : in Solution_List ) is
  begin
    put_line(file,p);
    Write_Solutions(file,sols);
  end put_line;

end QuadDobl_System_and_Solutions_io;
