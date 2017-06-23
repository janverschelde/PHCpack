with Communications_with_User;           use Communications_with_User;
with File_Scanning;                      use File_Scanning;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;

package body Prompt_for_Solutions is

  procedure Scan_Solutions
              ( file : in out file_type; onfile : in boolean;
                sols : in out Standard_Complex_Solutions.Solution_List;
                found : out boolean ) is
  begin
    if onfile then
      Scan_and_Skip(file,"THE SOLUTIONS",found);
      if found
       then get(file,sols);
       end if;
       Close(file);
    else
      found := false;
    end if;
  exception
    when others
      => put_line("Something is wrong with the solutions, will ignore...");
         Close(file);
         found := false;
  end Scan_Solutions;

  procedure Scan_Solutions
              ( file : in out file_type; onfile : in boolean;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                found : out boolean ) is
  begin
    if onfile then
      Scan_and_Skip(file,"THE SOLUTIONS",found);
      if found
       then get(file,sols);
       end if;
       Close(file);
    else
      found := false;
    end if;
  exception
    when others
      => put_line("Something is wrong with the solutions, will ignore...");
         Close(file);
         found := false;
  end Scan_Solutions;

  procedure Scan_Solutions
              ( file : in out file_type; onfile : in boolean;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                found : out boolean ) is
  begin
    if onfile then
      Scan_and_Skip(file,"THE SOLUTIONS",found);
      if found
       then get(file,sols);
       end if;
       Close(file);
    else
      found := false;
    end if;
  exception
    when others
      => put_line("Something is wrong with the solutions, will ignore...");
         Close(file);
         found := false;
  end Scan_Solutions;

  procedure Read_Solutions
              ( file : in out file_type; onfile : in boolean;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

    found : boolean;

  begin
    Scan_Solutions(file,onfile,sols,found);
    if not found then
      new_line;
      put_line("Reading the name of the file for the solutions.");
      Read_Name_and_Open_File(file);
      get(file,sols);
      Close(file);
    end if;
  end Read_Solutions;

  procedure Read_Solutions
              ( file : in out file_type; onfile : in boolean;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

    found : boolean;

  begin
    Scan_Solutions(file,onfile,sols,found);
    if not found then
      new_line;
      put_line("Reading the name of the file for the solutions.");
      Read_Name_and_Open_File(file);
      get(file,sols);
      Close(file);
    end if;
  end Read_Solutions;

  procedure Read_Solutions
              ( file : in out file_type; onfile : in boolean;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

    found : boolean;

  begin
    Scan_Solutions(file,onfile,sols,found);
    if not found then
      new_line;
      put_line("Reading the name of the file for the solutions.");
      Read_Name_and_Open_File(file);
      get(file,sols);
      Close(file);
    end if;
  end Read_Solutions;

end Prompt_for_Solutions;
