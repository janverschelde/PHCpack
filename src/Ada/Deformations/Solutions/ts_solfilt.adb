with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Drivers_for_Solution_Filters;      use Drivers_for_Solution_Filters;

procedure ts_solfilt is

-- DESCRIPTION :
--   Calls the driver to filter solutions.

  procedure Main is

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
    Driver_for_Solution_Filters(file,sols);
  end Main;

begin
  Main;
end ts_solfilt;
