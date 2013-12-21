with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;

procedure ts_commuser is

-- DESCRIPTION :
--   Tests the dialogues with the user.

  procedure Create_and_Write is

    file : file_type;
    c : character;

  begin
    put_line("Reading the name of an output file.");
    Read_Name_and_Create_File(file);
    put("Give a character : "); get(c);
    put(file,c);
    Close(file);
  end Create_and_Write;

  procedure Open_and_Read is

    file : file_type;
    c : character;

  begin
    put_line("Reading the name of an input file.");
    Read_Name_and_Open_File(file);
    get(file,c);
    put("The first character read : "); put(c); new_line;
    Close(file);
  end Open_and_Read;

  procedure Open_and_Append is

    file : file_type;
    c : character;

  begin
    put_line("Reading the name of a file to append to.");
    Read_Name_and_Append_File(file);
    put("Give a character : "); get(c);
    put(file,c);
    Close(file);
  end Open_and_Append;

  procedure Main is
  
    ans : character;

  begin
    new_line;
    put_line("Testing interactions with user to handle files.");
    new_line;
    loop
      put_line("Choose one of the following :");
      put_line("  0. Exit: leave this menu.");
      put_line("  1. Create an output file and write a character to it.");
      put_line("  2. Open a file for input and read a first character.");
      put_line("  3. Open a file and append a character.");
      put("Make your choice : "); Ask_Alternative(ans,"0123");
      new_line;
      case ans is
        when '1' => Create_and_Write;
        when '2' => Open_and_Read;
        when '3' => Open_and_Append;
        when others => exit;
      end case;
    end loop;
  end Main;

begin
  Main;
end ts_commuser;
