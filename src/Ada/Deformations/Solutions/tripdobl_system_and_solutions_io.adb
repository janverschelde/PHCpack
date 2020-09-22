with Communications_with_User;           use Communications_with_User;
with File_Scanning;                      use File_Scanning;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with TripDobl_Complex_Polynomials;
with TripDobl_Complex_Poly_Systems_io;   use TripDobl_Complex_Poly_Systems_io;
with TripDobl_Complex_Laurentials;
with TripDobl_Complex_Laur_Systems_io;   use TripDobl_Complex_Laur_Systems_io;
with TripDobl_Complex_Solutions_io;      use TripDobl_Complex_Solutions_io;

package body TripDobl_System_and_Solutions_io is

-- AUXILIARY ROUTINES :

  procedure Scan_for_Solutions
              ( file : in file_type; sols : out Solution_List;
                banner : in string := "SOLUTIONS" ) is

    found : boolean;

  begin
    Scan_and_Skip(file,banner,found);
    if found
     then get(file,sols);
    end if;
  end Scan_for_Solutions;

  procedure Write_Solutions
              ( file : in file_type; sols : in Solution_List;
                banner : in string := "THE SOLUTIONS :" ) is
  begin
    if not Is_Null(sols) then
      new_line(file);
      put_line(file,banner);
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    end if;
  end Write_Solutions;

-- TARGET ROUTINES :

  procedure get ( file : in file_type;
                  p : out Link_to_Poly_Sys; sols : out Solution_List;
                  banner : in string := "SOLUTIONS" ) is

  begin
    get(file,p);
    Scan_for_Solutions(file,sols,banner);
  end get;

  procedure get ( p : out Link_to_Poly_Sys; sols : out Solution_List;
                  banner : in string := "SOLUTIONS" ) is

    file : file_type;

  begin
    new_line;
    put_line("Reading the name of the input file ...");
    Read_Name_and_Open_File(file);
    get(file,p,sols,banner);
    close(file);
  end get;

  procedure get ( file : in file_type;
                  p : out Link_to_Laur_Sys; sols : out Solution_List;
                  banner : in string := "SOLUTIONS" ) is

  begin
    get(file,p);
    Scan_for_Solutions(file,sols,banner);
  end get;

  procedure get ( p : out Link_to_Laur_Sys; sols : out Solution_List;
                  banner : in string := "SOLUTIONS" ) is

    file : file_type;

  begin
    new_line;
    put_line("Reading the name of the input file ...");
    Read_Name_and_Open_File(file);
    get(file,p,sols,banner);
    close(file);
  end get;

  procedure put ( file : in file_type;
                  p : in Poly_Sys; sols : in Solution_List;
                  banner : in string := "THE SOLUTIONS :" ) is

    nbvar : constant natural32
          := TripDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first));
    nbequ : constant natural32 := natural32(p'last);

  begin
    if nbequ = nbvar then
      put(file,nbequ,1); new_line(file);
    else
      put(file,nbequ,1); put(file,"  ");
      put(file,nbvar,1); new_line(file);
    end if;
    put(file,p);
    Write_Solutions(file,sols,banner);
  end put;

  procedure put ( file : in file_type;
                  p : in Laur_Sys; sols : in Solution_List;
                  banner : in string := "THE SOLUTIONS :" ) is

    nbvar : constant natural32
          := TripDobl_Complex_Laurentials.Number_of_Unknowns(p(p'first));
    nbequ : constant natural32 := natural32(p'last);

  begin
    if nbequ = nbvar then
      put(file,nbequ,1); new_line(file);
    else
      put(file,nbequ,1); put(file,"  ");
      put(file,nbvar,1); new_line(file);
    end if;
    put(file,p);
    Write_Solutions(file,sols,banner);
  end put;

  procedure put_line ( file : in file_type;
                       p : in Poly_Sys; sols : in Solution_List;
                       banner : in string := "THE SOLUTIONS :" ) is
  begin
    put_line(file,p);
    Write_Solutions(file,sols,banner);
  end put_line;

  procedure put_line ( file : in file_type;
                       p : in Laur_Sys; sols : in Solution_List;
                       banner : in string := "THE SOLUTIONS :" ) is
  begin
    put_line(file,p);
    Write_Solutions(file,sols,banner);
  end put_line;

end TripDobl_System_and_Solutions_io;
