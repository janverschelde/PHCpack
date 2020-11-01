with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;
with Standard_Dictionary_Solutions_io;

package body Test_Dictionary_Solutions is

  procedure Standard_PHCpack_to_Dictionary is

    use Standard_Complex_Solutions;
    use Standard_Complex_Solutions_io;

    outfile : file_type;
    sols : Solution_List;

  begin
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(outfile);
    new_line; 
    Read(sols);
   -- put(standard_output,Length_Of(sols),Head_Of(sols).n,sols);
    new_line;
    put_line("See the output file for results...");
    new_line;
    Standard_Dictionary_Solutions_io.put(outfile,sols);
  end Standard_PHCpack_to_Dictionary;

  procedure Standard_Dictionary_to_PHCpack is

    use Standard_Complex_Solutions;
    use Standard_Complex_Solutions_io;

    infile,outfile : file_type;
    sols : Solution_List;

  begin
    new_line;
    put_line("Reading the name of the input file.");
    Read_Name_and_Open_File(infile);
    Standard_Dictionary_Solutions_io.get(infile,sols);
    close(infile);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(outfile);
    new_line;
    put_line("See the output file for results...");
    new_line;
    if not Is_Null(sols)
     then put(outfile,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    end if;
  end Standard_Dictionary_to_PHCpack;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to test format conversions of solution lists:");
    put_line("  1. from PHCpack to dictionary format, in standard floats;");
    put_line("  2. from dictionary to PHCpack format, in standard floats;");
    put("Type 1 or 2 to select the type of conversion : ");
    Ask_Alternative(ans,"12");
    case ans is
      when '1' => Standard_PHCpack_to_Dictionary;
      when '2' => Standard_Dictionary_to_PHCpack;
      when others => null;
    end case;
  end Main;

end Test_Dictionary_Solutions;
