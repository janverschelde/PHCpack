with text_io;                           use text_io;
with File_Management;

procedure ts_filman is

-- DESCRIPTION :
--   Test on management of file_type objects.

  c : character;

begin
  new_line;
  File_Management.Create_Output_File;
  put_line(File_Management.Link_to_Output.all,"hello world!");
  File_Management.Close_Output_File;
  File_Management.Open_Input_File;
  get(File_Management.Link_to_Input.all,c);
  put_line("The first character on the file : " & c);
end ts_filman;
