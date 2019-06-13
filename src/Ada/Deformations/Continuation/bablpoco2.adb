with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with DoblDobl_BlackBox_Continuations;    use DoblDobl_BlackBox_Continuations;

procedure bablpoco2 ( targetname,startname,outfilename : in string;
                      verbose : in integer32 := 0 ) is

  targetfile,startfile,outfile : file_type;
  poco : duration;

begin
  if verbose > 0 then
    put("At verbose level "); put(verbose,1);
    put_line(", in bablpoco2 ...");
  end if;
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
end bablpoco2;
