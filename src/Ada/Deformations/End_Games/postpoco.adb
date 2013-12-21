with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with valipoco;

procedure postpoco is

  pocofile,resultfile : file_type;

begin
  new_line;
  put_line("Reading name of the output file of poco.");
  Read_Name_and_Open_File(pocofile);
  new_line;
  put_line("Reading name of output file.");
  Read_Name_and_Create_File(resultfile);
  new_line;
  put_line("See output file for results...");
  new_line;
  valipoco(pocofile,resultfile);
end postpoco;
