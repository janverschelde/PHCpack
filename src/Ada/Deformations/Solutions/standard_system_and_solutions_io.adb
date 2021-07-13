with Communications_with_User;          use Communications_with_User;
with File_Scanning;                     use File_Scanning;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;

package body Standard_System_and_Solutions_io is

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

  begin
    get(file,integer(n),integer(m),p);
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
                  p : out Link_to_Poly_Sys; sols : out Solution_List;
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

  procedure get ( file : in file_type;
                  p : out Link_to_Laur_Sys; sols : out Solution_List;
                  banner : in string := "SOLUTIONS" ) is
  begin
    get(file,p);
    Scan_for_Solutions(file,sols,banner);
  end get;

  procedure put ( file : in file_type;
                  p : in Poly_Sys; sols : in Solution_List;
                  banner : in string := "THE SOLUTIONS :" ) is

    nbvar : constant natural32
          := Standard_Complex_Polynomials.Number_of_Unknowns(p(p'first));
    nbequ : constant natural32 := natural32(p'last);

  begin
    if nbequ = nbvar
     then put(file,nbequ,p);
     else put(file,nbequ,nbvar,p);
    end if;
    Write_Solutions(file,sols,banner);
  end put;

  procedure put ( file : in file_type;
                  p : in Laur_Sys; sols : in Solution_List;
                  banner : in string := "THE SOLUTIONS :" ) is

    nbvar : constant natural32
          := Standard_Complex_Laurentials.Number_of_Unknowns(p(p'first));
    nbequ : constant natural32 := natural32(p'last);

  begin
    if nbequ = nbvar
     then put(file,nbequ,p);
     else put(file,nbequ,nbvar,p);
    end if;
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

  procedure Scan_for_Start_System 
              ( infile : file_type; name : in Link_to_String;
                q : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                qsols : out Standard_Complex_Solutions.Solution_List;
                found : out boolean; verbose : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    len,dim : natural32 := 0;

  begin
    if vrblvl > 0 then
      put("-> in Standard_System_and_Solutions_io.");
      put_line("Scan_for_Start_System ...");
    end if;
    File_Scanning.Scan_and_Skip(infile,"START SYSTEM",found);
    if not found then
      if verbose then
        new_line;
        put_line("no start system found in " & name.all);
      end if;
    else
      get(infile,q);
      File_Scanning.Scan_and_skip(infile,"START SOLUTIONS",found);
      if not found then
        if verbose then
          new_line;
          put_line("no start solutions found in " & name.all);
        end if;
      else
        get(infile,len);
        get(infile,dim);
        get(infile,len,dim,qsols);
        if verbose then
          new_line;
          put("Read "); put(len,1);
          put(" solutions of dimension "); put(dim,1); put_line(".");
        end if;
      end if;
    end if;
  end Scan_for_Start_System;

  procedure Write_Scanned_Start_System
              ( name : in Link_to_String;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List ) is
  begin
    Write_Scanned_Start_System(standard_output,name,p,sols);
  end Write_Scanned_Start_System;

  procedure Write_Scanned_Start_System
              ( file : in file_type;
                name : in Link_to_String;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List ) is
  begin
    put(file,natural32(p'last),p);
    new_line(file);
    put_line(file,"TITLE : start system in file " & name.all);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
  end Write_Scanned_Start_System;

  procedure Main ( infilename,outfilename : in string;
                   vrblvl : in integer32 := 0 ) is

    infile,outfile : file_type;
    name : Link_to_String;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    found : boolean;
    verbose : constant boolean := (vrblvl > 0);

  begin
    if vrblvl > 0 then
      put_line("-> in standard_system_and_solutions_io.Main ...");
    end if;
    if infilename = "" then
      new_line;
      put_line("Reading the name of an input file ...");
      Read_Name_and_Open_File(infile,name);
    else
      Open_Input_File(infile,infilename,name);
    end if;
    Scan_for_Start_System(infile,name,lp,sols,found,verbose,vrblvl-1);
    if not found then
      if vrblvl > 0
       then put_line("no start system found in " & name.all);
      end if;
    else
      if outfilename = "" then
       -- new_line;
       -- put_line("Reading the name of an output file ...");
       -- Read_Name_and_Create_File(outfile);
        Write_Scanned_Start_System(name,lp.all,sols);
      else
        Create_Output_File(outfile,outfilename);
        Write_Scanned_Start_System(outfile,name,lp.all,sols);
        close(outfile);
      end if;
    end if;
  end Main;

end Standard_System_and_Solutions_io;
