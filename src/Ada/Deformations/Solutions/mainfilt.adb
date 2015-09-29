with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
--with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
--with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Drivers_for_Solution_Filters;      use Drivers_for_Solution_Filters;
with Standard_Select_Solutions;         use Standard_Select_Solutions;
-- with Drivers_for_Condition_Tables;      use Drivers_for_Condition_Tables;

procedure mainfilt ( infilename,outfilename : in string ) is

-- DESCRIPTION :
--   Calls the driver to filter solution lists.
--   There are two types of filters.  By default, expecting huge lists,
--   the scanning filter should be used, but for testing purposes, also
--   the "to list" filter can be enabled.

 -- procedure Read_Solutions ( sols : out Solution_List ) is
--
  -- DESCRIPTION :
  --   Prompts the user to provide a file name for the solutions,
  --   which are then loaded entirely into main memory.
  --   This may not work very well when the list is huge...
--
 --   infile : file_type;
 --   fail : boolean;
--
 -- begin
 --   loop
 --     new_line;
 --     put_line("Reading the name of the file for the solutions.");
 --     Read_Name_and_Open_File(infile);
 --     Prompt_to_Scan_Banner(infile,fail);
 --     exit when not fail;
 --     Close(infile);
 --     put("Something was wrong with the input file.");
 --     put_line("  Please try again...");
 --   end loop;
 --   get(infile,sols);
 -- exception
 --   when others => close(infile); raise;
 -- end Read_Solutions;

  procedure Read_Dimensions
              ( infile : out file_type; len,dim : out natural32 ) is

  -- DESCRIPTION :
  --   Prompts the user for the name of a file for the solutions,
  --   eventually preceeded by a polynomial system.  Only the dimensions 
  --   of the solution list are read, which is fine for huge lists.

    bannered,fail : boolean;

  begin
    loop
      new_line;
      put_line("Reading the name of the file for the solutions.");
      Read_Name_and_Open_File(infile);
      Prompt_to_Scan_Banner(infile,bannered,fail);
      exit when not fail;
      Close(infile);
      put("Something was wrong with the input file.");
      put_line("  Please try again...");
    end loop;
    Read_Dimensions(infile,len,dim,fail);
    put("Ready to process "); put(len,1);
    put(" solutions of dimension "); put(dim,1); put_line(" ... ");
  exception
    when others => close(infile); raise;
  end Read_Dimensions;

  function Prompt_for_Precision return integer32 is

  -- DESCRIPTION :
  --   Prompts the user to select the precision
  --   and returns 0 for double, 1 for double double,
  --   or 2 for quad double.

    ans : character;

  begin
    new_line;
    put_line("MENU for the precision : ");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => return 0;
      when '1' => return 1;
      when '2' => return 2;
      when others => return -1;
    end case;
  end Prompt_for_Precision;

  procedure Main is

    infile,outfile : file_type;
   -- list_filter : boolean := true;  -- use only for testing
   -- list_filter : boolean := false;   -- permits huge lists
   -- sols : Solution_List;
    len,dim : natural32;
    precision : integer32;

  begin
    new_line;
    put_line("Filtering solution lists subject to criteria.");
    declare
    begin
      --if list_filter
      -- then Read_Solutions(sols);
      -- else
      Read_Dimensions(infile,len,dim);
      --end if;
    exception
      when others
        => put_line("Please try again...");
          -- if list_filter
          --  then Read_Solutions(sols);
          --  else
           Read_Dimensions(infile,len,dim);
          -- end if;
    end;
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(outfile);
   -- if list_filter
   --  then Driver_for_Solution_Filters(outfile,sols);
   --  else
    precision := Prompt_for_Precision;
    case precision is
      when 0 => Driver_for_Standard_Solution_Filters(infile,outfile,len,dim);
      when 1 => Driver_for_DoblDobl_Solution_Filters(infile,outfile,len,dim);
      when 2 => Driver_for_QuadDobl_Solution_Filters(infile,outfile,len,dim);
      when others => null;
    end case;
   -- end if;
  end Main;

begin
  Main;
end mainfilt;
