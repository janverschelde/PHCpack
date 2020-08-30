with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Penta_Double_Vectors;
with Penta_Double_Vectors_io;
with PentDobl_Complex_Vectors;
with PentDobl_Complex_Vectors_io;
with PentDobl_Random_Vectors;

package body Test_PentDobl_Vectors is

  procedure Test_Penta_Double_Vectors is

    use Penta_Double_Vectors;
    use Penta_Double_Vectors_io;

    n : integer32 := 0;

  begin
    new_line;
    put("Give the dimension : "); get(n);
    declare
      v : constant Vector(1..n)
        := PentDobl_Random_Vectors.Random_Vector(1,n);
    begin
      put_line("a random vector : "); put_line(v);
    end;
  end Test_Penta_Double_Vectors;

  procedure Test_Penta_Double_Complex_Vectors is

    use PentDobl_Complex_Vectors;
    use PentDobl_Complex_Vectors_io;

    n : integer32 := 0;

  begin
    new_line;
    put("Give the dimension : "); get(n);
    declare
      v : constant Vector(1..n)
        := PentDobl_Random_Vectors.Random_Vector(1,n);
    begin
      put_line("a random vector : "); put_line(v);
    end;
  end Test_Penta_Double_Complex_Vectors;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Interactive testing of vectors of penta doubles.");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line("  0. exit this program;");
      put_line("  1. random vector of penta double numbers;");
      put_line("  2. random vector of complex penta double numbers.");
      put("Type 1 or 2 to select a test, or 0 to exit : ");
      Ask_Alternative(ans,"012");
      exit when (ans = '0');
      case ans is
        when '1' => Test_Penta_Double_Vectors;
        when '2' => Test_Penta_Double_Complex_Vectors;
        when others => null;
      end case;
    end loop;
  end Main;

end Test_PentDobl_Vectors;
