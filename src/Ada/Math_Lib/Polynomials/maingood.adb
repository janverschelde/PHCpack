with text_io;                            use text_io;
with Parse_Strings_to_Polynomials;       use Parse_Strings_to_Polynomials;

procedure maingood ( infilename,outfilename : in string ) is

  procedure Main is
  begin
    if infilename /= "" then
      Read_from_File(infilename,outfilename);
    else
      new_line;
      put_line("Reading the name of the input file ...");
      Read_Input_File_Name; 
    end if;
  end Main;

begin
  Main;
end maingood;
