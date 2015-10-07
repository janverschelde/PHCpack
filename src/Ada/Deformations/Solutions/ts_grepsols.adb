with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors;           use Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_Solution_Filters;          use Standard_Solution_Filters;
with Standard_Select_Solutions;          use Standard_Select_Solutions;

procedure ts_grepsols is

-- DESCRIPTION :
--   Reads a solution list from file.
--   Writes a selection of a solution list to file.

  function Read_Selection return Vector is

  -- DESCRIPTION :
  --   Returns a vector of natural numbers.

    n : integer32 := 0;

  begin
    put("Give the number of solutions to select : "); get(n);
    declare
      v : Vector(1..n);
    begin
      put("Give "); put(n,1); put(" natural increasing numbers : ");
      get(v);
      return v;
    end;
  end Read_Selection;

  procedure Select_from_List ( infile,outfile : in file_type ) is

  -- DESCRIPTION :
  --   Does a selection of the solutions by reading in the whole list.
  --   This works fine for relatively small lists of solutions.

    sols : Solution_List;
    bannered,fail : boolean;

  begin
    Prompt_to_Scan_Banner(infile,bannered,fail);
    get(infile,sols);
    new_line;
    put("There are "); put(Length_Of(sols),1); put_line(" solutions given.");
    declare
      sel : constant Vector := Read_Selection;
      selsols : Solution_List := Select_Solutions(sols,sel);
    begin
      put_line("The selected solutions : "); put(selsols);
      put(outfile,Length_Of(selsols),natural32(Head_Of(selsols).n),selsols);
    end;
  end Select_from_List;

  procedure Select_by_Scanning ( infile,outfile : in file_type ) is

  -- DESCRIPTION :
  --   Scans the solution lists one after the other and does not
  --   store all solutions in main memory.

    bannered,fail : boolean;
    len,dim,cnt : natural32;

  begin
    Prompt_to_Scan_Banner(infile,bannered,fail);
    Read_Dimensions(infile,len,dim,fail);
    put("Ready to process "); put(len,1);
    put(" solutions of dimension "); put(dim,1); put_line(" ... ");
    declare
      sel : constant Vector := Read_Selection;
    begin
      Select_Solutions(infile,outfile,len,dim,sel,cnt);
    end;
    put("Wrote "); put(cnt,1); put_line(" solutions to file.");
  end Select_by_Scanning;

  procedure Main is

    infile,outfile : file_type;

  begin
    new_line;
    put_line("Selecting elements from a list of solutions.");
    new_line;
    put_line("Reading the name of the file with the solutions.");
    Read_Name_and_Open_File(infile);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(outfile);
   -- Select_from_List(infile,outfile);
    Select_by_Scanning(infile,outfile);
    Close(infile); Close(outfile);
  end Main;

begin
  Main;
end ts_grepsols;
