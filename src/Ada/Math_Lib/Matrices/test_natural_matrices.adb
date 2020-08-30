with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Natural_Matrices;
with Standard_Natural_Matrices_io;
with Multprec_Natural_Matrices;
with Multprec_Natural_Matrices_io;

package body Test_Natural_Matrices is

  procedure Test_Standard_io is

    use Standard_Natural_Matrices;
    use Standard_Natural_Matrices_io;

    n,m : integer32 := 0;

  begin
    put("Give the number of rows : "); get(n);
    put("Give the number of columns : "); get(m);
    declare
      mat : Matrix(1..n,1..m);
    begin
      put("Give "); put(n,1); put("x"); put(m,1);
      put_line(" natural matrix : "); get(mat);
      put_line("Your matrix : "); put(mat); new_line;
    end;
  end Test_Standard_io;

  procedure Test_Multprec_io is

    use Multprec_Natural_Matrices;
    use Multprec_Natural_Matrices_io;

    n,m : integer32 := 0;

  begin
    put("Give the number of rows : "); get(n);
    put("Give the number of columns : "); get(m);
    declare
      mat : Matrix(1..n,1..m);
    begin
      put("Give "); put(n,1); put("x"); put(m,1);
      put_line(" natural matrix : "); get(mat);
      put_line("Your matrix : "); put(mat); new_line;
    end;
  end Test_Multprec_io;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Interactive testing of matrices of natural numbers");
    new_line;
    loop
      put_line("Choose one of the following : ");
      put_line("  1. io of matrices of standard natural numbers.");
      put_line("  2. io of matrices of multi-precision natural numbers.");
      put("Make your choice (1/2) : "); get(ans);
      case ans is
        when '1' => Test_Standard_io;
        when '2' => Test_Multprec_io;
        when others => null;
      end case;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Main;

end Test_Natural_Matrices;
