with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_Vectors;           use Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Complex_Solutions;
with Standard_Select_Solutions;          use Standard_Select_Solutions;
with DoblDobl_Complex_Solutions;
with DoblDobl_Select_Solutions;          use DoblDobl_Select_Solutions;
with QuadDobl_Complex_Solutions;
with QuadDobl_Select_Solutions;          use QuadDobl_Select_Solutions;

procedure ts_selsol is

-- DESCRIPTION :
--   Selecting a list of solutions from file.

  procedure Read_File_and_Scan_Solutions 
              ( sv : in Vector; dim : out natural32;
                sa : out Standard_Complex_Solutions.Solution_Array ) is

  -- DESCRIPTION :
  --   Prompts the user for a file name for a solution list
  --   and scans the file for a selection of the solutions in the list.
  --   The solutions on return are in standard double precision.

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
                sa : out DoblDobl_Complex_Solutions.Solution_Array ) is

  -- DESCRIPTION :
  --   Prompts the user for a file name for a solution list
  --   and scans the file for a selection of the solutions in the list.
  --   The solutions on return are in double double precision.

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
                sa : out QuadDobl_Complex_Solutions.Solution_Array ) is

  -- DESCRIPTION :
  --   Prompts the user for a file name for a solution list
  --   and scans the file for a selection of the solutions in the list.
  --   The solutions on return are in quad double precision.

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
                sols : out Standard_Complex_Solutions.Solution_List ) is

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

  procedure Standard_Scan_for_Selection
              ( n : in natural32; rv : in Vector ) is

  -- DESCRIPTION :
  --   Sorts the sequence of n numbers in the vector rv
  --   and prompts the user for a file name.

    sv : constant Vector(rv'range) := Standard_Select_Solutions.Sort(rv);
    sa : Standard_Complex_Solutions.Solution_Array(rv'range);
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
  end Standard_Scan_for_Selection;

  procedure DoblDobl_Scan_for_Selection
              ( n : in natural32; rv : in Vector ) is

  -- DESCRIPTION :
  --   Sorts the sequence of n numbers in the vector rv
  --   and prompts the user for a file name.

    sv : constant Vector(rv'range) := DoblDobl_Select_Solutions.Sort(rv);
    sa : DoblDobl_Complex_Solutions.Solution_Array(rv'range);
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
  end DoblDobl_Scan_for_Selection;

  procedure QuadDobl_Scan_for_Selection
              ( n : in natural32; rv : in Vector ) is

  -- DESCRIPTION :
  --   Sorts the sequence of n numbers in the vector rv
  --   and prompts the user for a file name.

    sv : constant Vector(rv'range) := QuadDobl_Select_Solutions.Sort(rv);
    sa : QuadDobl_Complex_Solutions.Solution_Array(rv'range);
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
  end QuadDobl_Scan_for_Selection;

  function Get_Selection return Link_to_Vector is

  -- DESCRIPTION :
  --   Prompts the user for a select of indices,
  --   either from file or interactively typed in.
  --
    res : Link_to_Vector;
    num_file : file_type;
    n : natural32 := 0;
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
      get(num_file,n,res);
    else
      new_line;
      put("Give the number of indices : "); get(n);
      res := new Vector(1..integer32(n));
      put("Give "); put(n,1); put(" indices : ");
      for i in 1..integer32(n) loop
        get(res(i));
      end loop;
      skip_line; -- skip last newline symbol
      put("-> your indices : "); put(res); new_line;
    end if;
    return res;
  end Get_Selection;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a selection of indices
  --   and the precision before getting solutions from file.

    rv : Link_to_Vector := Get_Selection;
    nb : constant natural32 := natural32(rv'last);
    ans : character;

  begin
    new_line;
    put_line("MENU to select the precision : ");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision;");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Scan_for_Selection(nb,rv.all);
      when '1' => DoblDobl_Scan_for_Selection(nb,rv.all);
      when '2' => QuadDobl_Scan_for_Selection(nb,rv.all);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_selsol;
