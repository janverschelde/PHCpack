with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with File_Scanning;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;


procedure ts_trasols_io is

-- DESCRIPTION :
--   Development of processing of the output file of a path tracker.

  procedure Read ( file : in file_type; verbose : in boolean := false ) is

  -- DESCRIPTION :
  --   Reads the target system, the start system, the start solution,
  --   and the solutions of the target system (in this order) from
  --   the file.  If verbose, then extra output is written to screen.

    lp,lq : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    found : boolean;
    qsols,psols : Standard_Complex_Solutions.Solution_List;

  begin
    get(file,lp);
    if verbose then
      new_line;
      put_line("The target system :");
      put(lp.all);
    end if;
    File_Scanning.Scan_and_Skip(file,"START SYSTEM",found);
    if found then
      get(file,lq);
      if verbose then
        new_line;
        put_line("The start system :");
        put(lq.all);
      end if;
      File_Scanning.Scan_and_Skip(file,"START SOLUTIONS",found);
      if found then
        get(file,qsols);
        if verbose then
          new_line;
          put("Read ");
          put(Standard_Complex_Solutions.Length_Of(qsols),1);
          put_line(" start solutions.");
        end if;
        File_Scanning.Scan_and_Skip(file,"SOLUTIONS",found);
        if found then
          get(file,psols);
          if verbose then
            new_line;
            put("Read ");
            put(Standard_Complex_Solutions.Length_Of(psols),1);
            put_line(" solutions.");
          end if;
        end if;
      end if;
    end if;
  end Read;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a file name, the output of a path tracker.

    infile : file_type;

  begin
    new_line;
    put_line("Reading the name of a file with output of a tracker ...");
    Read_Name_and_Open_File(infile);
    new_line;
    put_line("Reading the contents of the file ...");
    Read(infile,true);
  end Main;

begin
  Main;
end ts_trasols_io;
