with text_io;                           use text_io;
with Communications_with_User;
with demics_input_data;

package body Test_DEMiCs_Input is

  procedure read_data_from_file is
  begin
    null;
  end read_data_from_file;

  procedure interactive_input_data is
  begin
    null;
  end interactive_input_data;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Testing input data for DEMiCs ...");
    new_line;
    put("Reading data from file ? (y/n) ");
    Communications_with_User.Ask_Yes_or_No(ans);
    if ans = 'y'
     then read_data_from_file;
     else interactive_input_data;
    end if;
  end Main;

end Test_DEMiCs_Input;
