with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Double_Double_Matrices;
with Double_Double_Matrices_io;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Matrices_io;
with DoblDobl_Random_Matrices;

package body Test_DoblDobl_Matrices is

  procedure Test_Double_Double_Matrices is

    use Double_Double_Matrices;
    use Double_Double_Matrices_io;

    n,m : natural32 := 0;

  begin
    put("Give the number of rows : "); get(n);
    put("Give the number of columns : "); get(m);
    declare
      A : constant Matrix(1..integer32(n),1..integer32(m)) 
        := DoblDobl_Random_Matrices.Random_Matrix(n,m);
    begin
      put_line("a random matrix : "); put(A);
    end;
  end Test_Double_Double_Matrices;

  procedure Test_Double_Double_Complex_Matrices is

  -- DESCRIPTION :
  --   Prompts the user for number of rows and columns
  --   and then shows a matrix of random complex double doubles.

    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Matrices_io;

    n,m : natural32 := 0;

  begin
    put("Give the number of rows : "); get(n);
    put("Give the number of columns : "); get(m);
    declare
      A : constant Matrix(1..integer32(n),1..integer32(m))
        := DoblDobl_Random_Matrices.Random_Matrix(n,m);
    begin
      put_line("a random matrix : "); put(A);
    end;
  end Test_Double_Double_Complex_Matrices;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Interactive testing of matrices of double doubles.");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line("  0. exit this program;");
      put_line("  1. random matrix of double double numbers;");
      put_line("  2. random matrix of complex double double numbers;");
      put("Type 0, 1, or 2 to make your choice : ");
      Ask_Alternative(ans,"012");
      exit when (ans = '0');
      case ans is
        when '1' => Test_Double_Double_Matrices;
        when '2' => Test_Double_Double_Complex_Matrices;
        when others => null;
      end case;
    end loop;
  end Main;

end Test_DoblDobl_Matrices;
