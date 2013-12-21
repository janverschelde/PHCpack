with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_Vectors;           use Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Drivers_for_Condition_Tables;       use Drivers_for_Condition_Tables;
with Standard_Select_Solutions;          use Standard_Select_Solutions;

procedure ts_selsol is

-- DESCRIPTION :
--   Selecting a list of solutions from file.

  procedure Read_File_and_Scan_Solutions 
              ( sv : in Vector; dim : out natural32;
                sa : out Solution_Array ) is

  -- DESCRIPTION :
  --   Prompts the user for a file name for a solution list
  --   and scans the file for a selection of the solutions in the list.

  -- ON ENTRY :
  --   sv       selection of solutions, sorted in increasing order.

  -- ON RETURN :
  --   dim      length of the solution vectors;
  --   sa       selected array of solutions, as they are read from file,
  --            according to the numbers in sv.

    file : file_type;
    length : natural32;
    bannered,fail : boolean;

  begin
    new_line;
    put_line("Reading the name of the input file for the solutions...");
    Read_Name_and_Open_File(file);
    Scan_Banner_Dimensions(file,length,dim,bannered,fail);
    if fail then
      put("Format of the solution list on file is incorrect.");
      put_line("  Please try again.");
    else
      put("Ready to scan "); put(length,1);
      put(" solutions of dimension "); put(dim,1); put_line(".");
      new_line;
      Scan_Solutions(file,length,dim,sv,sa);
    end if;
  end Read_File_and_Scan_Solutions;

  procedure Read_File_and_Scan_Solutions 
              ( sv : in Vector; dim : out natural32;
                sols : out Solution_List ) is

  -- DESCRIPTION :
  --   Prompts the user for a file name for a solution list
  --   and scans the file for a selection of the solutions in the list.

  -- ON ENTRY :
  --   sv       selection of solutions, sorted in increasing order.

  -- ON RETURN :
  --   dim      length of the solution vectors;
  --   sols     selected solutions, as they are read from file,
  --            according to the numbers in sv.

    file : file_type;
    length : natural32;
    bannered,fail : boolean;

  begin
    new_line;
    put_line("Reading the name of the input file for the solutions...");
    Read_Name_and_Open_File(file);
    Scan_Banner_Dimensions(file,length,dim,bannered,fail);
    if fail then
      put("Format of the solution list on file is incorrect.");
      put_line("  Please try again.");
    else
      put("Ready to scan "); put(length,1);
      put(" solutions of dimension "); put(dim,1); put_line(".");
      new_line;
      Scan_Solutions(file,length,dim,sv,sols);
    end if;
  end Read_File_and_Scan_Solutions;

  procedure Scan_for_Selection ( n : in natural32; rv : in Vector ) is

  -- DESCRIPTION :
  --   Sorts the sequence of n numbers in the vector rv
  --   and prompts the user for a file name.

    sv : constant Vector(rv'range) := Sort(rv);
    sa : Solution_Array(rv'range);
    dim : natural32;
    file : file_type; 

  begin
    new_line;
    put("The number of solutions to be read : "); put(n,1); new_line;
    put("Numbers to solutions : "); put(rv); new_line;
    put("The sorted numbers : "); put(sv); new_line;
    Read_File_and_Scan_Solutions(sv,dim,sa);
    new_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(file);
    Write_Selection(file,dim,rv,sv,sa);
  end Scan_for_Selection;

  procedure Main is

    num_file : file_type;
    n : natural32 := 0;
    rv : Link_to_Vector;
    ans : character;

  begin
    new_line;
    put_line("Selecting a list of solutions from file.");
    new_line;
    put_line("MENU to read the index numbers of the solutions : ");
    put_line("  1. provide a file name with numbers; or");
    put_line("  2. interactively type in the numbers.");
    put("Type 1 or 2 to make a choice : ");
    Ask_Alternative(ans,"12");
    if ans = '1' then
      new_line;
      put_line("Reading the name of the file for the numbers...");
      Read_Name_and_Open_File(num_file);
      get(num_file,n);
      get(num_file,n,rv);
    else
      new_line;
      put("Give the number of indices : "); get(n);
      rv := new Vector(1..integer32(n));
      put("Give "); put(n,1); put(" indices : ");
      for i in 1..integer32(n) loop
        get(rv(i));
      end loop;
      put("-> your indices : "); put(rv); new_line;
    end if;
    Scan_for_Selection(n,rv.all);
  end Main;

begin
  Main;
end ts_selsol;
