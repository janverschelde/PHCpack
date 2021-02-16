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

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the name of an output file
  --   and scans for the start system and start solutions.

    infile,outfile : file_type;
    found : boolean;
    name : Link_to_String;
    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    len,dim : natural32 := 0;

  begin
    new_line;
    put_line("Reading the name of an input file ...");
    Read_Name_and_Open_File(infile,name);
    File_Scanning.Scan_and_Skip(infile,"START SYSTEM",found);
    if not found then
      put_line("no start system found in " & name.all);
    else
      new_line;
      put_line("Reading the name of an output file ...");
      Read_Name_and_Create_File(outfile);
      get(infile,lp);
      put(outfile,natural32(lp'last),lp.all);
      new_line(outfile);
      put_line(outfile,"TITLE : start system in file " & name.all);
      File_Scanning.Scan_and_skip(infile,"START SOLUTIONS",found);
      if not found then
        put_line("no start solutions found in " & name.all);
      else
        get(infile,len);
        get(infile,dim);
        get(infile,len,dim,sols);
        put("Read "); put(len,1);
        put(" solutions of dimension "); put(dim,1); put_line(".");
        new_line(outfile);
        put_line(outfile,"THE SOLUTIONS :");
        put(outfile,len,dim,sols);
      end if;
      close(outfile);
    end if;
  end Main;

begin
  Main;
end ts_getstart;
