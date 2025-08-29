with text_io;                           use text_io;
with Communications_with_User;
with demics_input_data;
with demics_input_main;

package body Test_DEMiCs_Input is

  procedure read_data_from_file is

    use demics_input_data.class_dataSet;

    ans : character;
    data : dataSet;
    fail : boolean;

  begin
    put("Intermediate output wanted ? (y/n) ");
    Communications_with_User.Ask_Yes_or_No(ans);
    new_line;
    if ans = 'y'
     then demics_input_main.read_data_from_file(data,fail,1);
     else demics_input_main.read_data_from_file(data,fail,0);
    end if;
  end read_data_from_file;

  procedure interactive_input_data is

    use demics_input_data.class_dataSet;

    ans : character;
    data : dataSet;
    fail : boolean;

  begin
    put("Intermediate output wanted ? (y/n) ");
    Communications_with_User.Ask_Yes_or_No(ans);
    new_line;
    if ans = 'y'
     then demics_input_main.interactive_input_data(data,fail,1);
     else demics_input_main.interactive_input_data(data,fail,0);
    end if;
  end interactive_input_data;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Testing input data for DEMiCs ...");
    new_line;
    put("Reading data from file ? (y/n) ");
    Communications_with_User.Ask_Yes_or_No(ans);
    new_line;
    if ans = 'y'
     then read_data_from_file;
     else interactive_input_data;
    end if;
  end Main;

end Test_DEMiCs_Input;
