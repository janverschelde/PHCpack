with text_io;                            use text_io;
with String_Splitters;                   use String_Splitters;
with Communications_with_User;           use Communications_with_User;
with File_Scanning;                      use File_Scanning;

procedure ts_filescan is

-- DESCRIPTION :
--   Reads a file name and a banner, and scans for it.

  procedure Main is

    file : file_type;
    ans : character;
 
  begin
    new_line;
    put_line("Interactive test of file scanning.");
    new_line;
    loop
      put_line("Reading the name of the input file.");
      Read_Name_and_Open_File(file);
      loop
        put_line("Reading a banner.");
        declare
          banner : constant String := Read_String;
          found : boolean := false;
        begin
          Scan(file,banner,found);
          if found
           then put_line("The banner has been found.");
           else put_line("The banner has not been found.");
          end if; 
        end;
        put("Do you want to test more banners ? (y/n) ");
        Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
        Reset(file);
      end loop;
      Close(file);
      put("Do you want to scan other files ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Main;

begin
  Main;
end ts_filescan;
