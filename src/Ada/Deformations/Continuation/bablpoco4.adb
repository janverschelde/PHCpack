with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with QuadDobl_BlackBox_Continuations;    use QuadDobl_BlackBox_Continuations;

procedure bablpoco4 ( targetname,startname,outfilename : in string ) is

  targetfile,startfile,outfile : file_type;
  poco : duration;

begin
  if targetname /= "" then
    Open_Input_File(targetfile,targetname);
  else
    new_line;
    put_line("Reading the name of the file for the target system.");
    Read_Name_and_Open_File(targetfile);
  end if;
  if startname /= "" then
    Open_Input_File(startfile,startname);
  else
    new_line;
    put_line("Reading the name of the file for the start system.");
    Read_Name_and_Open_File(startfile);
  end if;
  Create_Output_File(outfile,outfilename);
  Black_Box_Polynomial_Continuation(targetfile,startfile,outfile,poco);
  Close(targetfile); Close(startfile); Close(outfile);
end bablpoco4;
