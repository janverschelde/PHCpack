with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with String_Splitters;                   use String_Splitters;
with File_Scanning;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;

procedure ts_getstart is

-- DESCRIPTION :
--   Tests the scanning of start system and start solutions
--   from an output file of the blackbox solver of phc.

  procedure Scan_for_Start_System 
              ( infile : file_type; name : in Link_to_String;
                q : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                qsols : out Standard_Complex_Solutions.Solution_List;
                found : out boolean; verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Scans the infile for a start system and start solutions.

  -- ON ENTRY :
  --   infile   must be opened for input;
  --   name     name of the infile for error messages;
  --   verbose  if true, then error messages are displayed
  --            when a start system or its solutions are not found.

  -- ON RETURN :
  --   q        a start system, if found;
  --   qsols    start solutions, if found;
  --   found    if true, then both q and qsols are found,
  --            otherwise, either q and/or qsols were not present.

    len,dim : natural32 := 0;

  begin
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

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the name of an output file
  --   and scans for the start system and start solutions.

    infile,outfile : file_type;
    name : Link_to_String;
    found : boolean;
    lp : Link_to_Poly_Sys;
    sols : Solution_List;

  begin
    new_line;
    put_line("Reading the name of an input file ...");
    Read_Name_and_Open_File(infile,name);
    Scan_for_Start_System(infile,name,lp,sols,found);
    if not found then
      put_line("no start system found in " & name.all);
    else
      new_line;
      put_line("Reading the name of an output file ...");
      Read_Name_and_Create_File(outfile);
      put(outfile,natural32(lp'last),lp.all);
      new_line(outfile);
      put_line(outfile,"TITLE : start system in file " & name.all);
      new_line(outfile);
      put_line(outfile,"THE SOLUTIONS :");
      put(outfile,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      close(outfile);
    end if;
  end Main;

begin
  Main;
end ts_getstart;
