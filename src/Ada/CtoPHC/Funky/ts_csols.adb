with text_io;                          use text_io;
with Communications_with_User;         use Communications_with_User;
with Standard_Natural_Numbers;         use Standard_Natural_Numbers;
with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with C_Integer_Arrays;                 use C_Integer_Arrays;
with C_Double_Arrays;                  use C_Double_Arrays;
with Standard_Complex_Solutions;       use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;    use Standard_Complex_Solutions_io;
with Coefficient_Solution_Vectors;     use Coefficient_Solution_Vectors;

procedure ts_csols is

-- DESCRIPTION :
--   Test on the coefficient representation of solutions.

  procedure Test_Conversions
              ( file : in file_type; sols : in Solution_List ) is

    n : constant integer32 := Head_Of(sols).n;
    m : constant C_Integer_Array := Multiplicities(sols);
    c : constant C_Double_Array := Coefficients(sols);
    s : Solution_List;

  begin
   -- put(file,"The Multiplicities : "); put(file,m); new_line(file);
   -- put_line(file,"The Coefficients : "); put_line(file,c);
    s := Create(natural32(n),m,c);
    put_line(file,"The created solution list :");
    put(file,Length_Of(s),natural32(Head_Of(s).n),s);
  end Test_Conversions;

  procedure Main is

    file : file_type;
    sols : Solution_List;

  begin
    new_line;
    put_line("Testing coefficient representations of solutions...");
    new_line;
    Read(sols);
    new_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(file);
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    new_line(file);
    Test_Conversions(file,sols);
  end Main;

begin
  Main;
end ts_csols;
