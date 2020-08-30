with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Octo_Double_Matrices;
with Octo_Double_Matrices_io;
with OctoDobl_Complex_Matrices;
with OctoDobl_Complex_Matrices_io;
with OctoDobl_Random_Matrices;

package body Test_OctoDobl_Matrices is

  procedure Test_Octo_Double_Matrices is

    use Octo_Double_Matrices;
    use Octo_Double_Matrices_io;

    n,m : natural32 := 0;

  begin
    put("Give the number of rows : "); get(n);
    put("Give the number of columns : "); get(m);
    declare
      A : constant Matrix(1..integer32(n),1..integer32(m)) 
        := OctoDobl_Random_Matrices.Random_Matrix(n,m);
    begin
      put_line("a random matrix : "); put(A);
    end;
  end Test_Octo_Double_Matrices;

  procedure Test_Octo_Double_Complex_Matrices is

    use OctoDobl_Complex_Matrices;
    use OctoDobl_Complex_Matrices_io;

    n,m : natural32 := 0;

  begin
    put("Give the number of rows : "); get(n);
    put("Give the number of columns : "); get(m);
    declare
      A : constant Matrix(1..integer32(n),1..integer32(m))
        := OctoDobl_Random_Matrices.Random_Matrix(n,m);
    begin
      put_line("a random matrix : "); put(A);
    end;
  end Test_Octo_Double_Complex_Matrices;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Interactive testing of matrices of octo doubles.");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line("  0. exit this program;");
      put_line("  1. random matrix of octo double numbers;");
      put_line("  2. random matrix of complex octo double numbers.");
      put("Type 0, 1, or 2 to make your choice : ");
      Ask_Alternative(ans,"012");
      exit when (ans = '0');
      case ans is
        when '1' => Test_Octo_Double_Matrices;
        when '2' => Test_Octo_Double_Complex_Matrices;
        when others => null;
      end case;
    end loop;
  end Main;

end Test_OctoDobl_Matrices;
