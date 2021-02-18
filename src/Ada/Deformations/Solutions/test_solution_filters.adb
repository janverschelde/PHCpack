with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;     use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;     use QuadDobl_Complex_Solutions_io;
with Main_Solution_Filters;             use Main_Solution_Filters;

package body Test_Solution_Filters is

  procedure Standard_Filter is

    use Standard_Complex_Solutions;

    file : file_type;
    sols : Solution_List;

  begin
    new_line;
    put_line("Filtering solution lists subject to criteria.");
    new_line;
    Read(sols);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    Main_Solution_Filter(file,sols);
  end Standard_Filter;

  procedure DoblDobl_Filter is

    use DoblDobl_Complex_Solutions;

    file : file_type;
    sols : Solution_List;

  begin
    new_line;
    put_line("Filtering solution lists subject to criteria.");
    new_line;
    Read(sols);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    Main_Solution_Filter(file,sols);
  end DoblDobl_Filter;

  procedure QuadDobl_Filter is

    use QuadDobl_Complex_Solutions;

    file : file_type;
    sols : Solution_List;

  begin
    new_line;
    put_line("Filtering solution lists subject to criteria.");
    new_line;
    Read(sols);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    Main_Solution_Filter(file,sols);
  end QuadDobl_Filter;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Filtering solution lists subject to criteria.");
    new_line;
    put_line("MENU for the precision : ");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Filter;
      when '1' => DoblDobl_Filter;
      when '2' => QuadDobl_Filter;
      when others => null;
    end case;
  end Main;

end Test_Solution_Filters;
