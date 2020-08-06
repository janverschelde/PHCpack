with Communications_with_User;          use Communications_with_User;
with File_Scanning;                     use File_Scanning;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Poly_Systems_io;  use Multprec_Complex_Poly_Systems_io;
with Multprec_Complex_Laurentials;
with Multprec_Complex_Laur_Systems_io;  use Multprec_Complex_Laur_Systems_io;
with Multprec_Complex_Solutions_io;     use Multprec_Complex_Solutions_io;

package body Multprec_System_and_Solutions_io is

  procedure get ( n,m : out natural32; p : out Link_to_Array_of_Strings;
                  sols : out Solution_List;
                  banner : in string := "SOLUTIONS" ) is

    file : file_type;

  begin
    new_line;
    put_line("Reading the name of the input file ...");
    Read_Name_and_Open_File(file);
    get(file,n,m,p,sols,banner);
    close(file);
  end get;

  procedure get ( file : in file_type;
                  n,m : out natural32; p : out Link_to_Array_of_Strings;
                  sols : out Solution_List;
                  banner : in string := "SOLUTIONS" ) is

    found : boolean;

  begin
    get(file,integer(n),integer(m),p);
    Scan_and_Skip(file,banner,found);
    if found
     then get(file,sols);
    end if;
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
                  p : out Link_to_Poly_Sys; sols : out Solution_List;
                  banner : in string := "SOLUTIONS" ) is

    found : boolean;

  begin
    get(file,p);
    Scan_and_Skip(file,banner,found);
    if found
     then get(file,sols);
    end if;
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

  procedure get ( file : in file_type;
                  p : out Link_to_Laur_Sys; sols : out Solution_List;
                  banner : in string := "SOLUTIONS" ) is

    found : boolean;

  begin
    get(file,p);
    Scan_and_Skip(file,banner,found);
    if found
     then get(file,sols);
    end if;
  end get;

  procedure put ( file : in file_type;
                  p : in Poly_Sys; sols : in Solution_List;
                  banner : in string := "THE SOLUTIONS :" ) is

    nbvar : constant natural32
          := Multprec_Complex_Polynomials.Number_of_Unknowns(p(p'first));
    nbequ : constant natural32 := natural32(p'last);

  begin
    if nbequ = nbvar then
      put(file,nbequ,1); new_line(file);
    else
      put(file,nbequ,1); put(file,"  ");
      put(file,nbvar,1); new_line(file);
    end if;
    put(file,p);
    if not Is_Null(sols) then
      new_line(file);
      put_line(file,banner);
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    end if;
  end put;

  procedure put ( file : in file_type;
                  p : in Laur_Sys; sols : in Solution_List;
                  banner : in string := "THE SOLUTIONS :" ) is

    nbvar : constant natural32
          := Multprec_Complex_Laurentials.Number_of_Unknowns(p(p'first));
    nbequ : constant natural32 := natural32(p'last);

  begin
    if nbequ = nbvar then
      put(file,nbequ,1); new_line(file);
    else
      put(file,nbequ,1); put(file,"  ");
      put(file,nbvar,1); new_line(file);
    end if;
    put(file,p);
    if not Is_Null(sols) then
      new_line(file);
      put_line(file,banner);
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    end if;
  end put;

  procedure put_line ( file : in file_type;
                       p : in Poly_Sys; sols : in Solution_List;
                       banner : in string := "THE SOLUTIONS :" ) is

    nbvar : constant natural32
          := Multprec_Complex_Polynomials.Number_of_Unknowns(p(p'first));
    nbequ : constant natural32 := natural32(p'last);

  begin
    if nbequ = nbvar then
      put(file,nbequ,1); new_line(file);
    else
      put(file,nbequ,1); put(file,"  ");
      put(file,nbvar,1); new_line(file);
    end if;
    put_line(file,p);
    if not Is_Null(sols) then
      new_line(file);
      put_line(file,banner);
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    end if;
  end put_line;

  procedure put_line ( file : in file_type;
                       p : in Laur_Sys; sols : in Solution_List;
                       banner : in string := "THE SOLUTIONS :" ) is

    nbvar : constant natural32
          := Multprec_Complex_Laurentials.Number_of_Unknowns(p(p'first));
    nbequ : constant natural32 := natural32(p'last);

  begin
    if nbequ = nbvar then
      put(file,nbequ,1); new_line(file);
    else
      put(file,nbequ,1); put(file,"  ");
      put(file,nbvar,1); new_line(file);
    end if;
    put_line(file,p);
    if not Is_Null(sols) then
      new_line(file);
      put_line(file,banner);
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    end if;
  end put_line;

end Multprec_System_and_Solutions_io;
