with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Multprec_Natural_Numbers;          use Multprec_Natural_Numbers;
with Multprec_Natural_Numbers_io;       use Multprec_Natural_Numbers_io;
with Multprec_Integer_Numbers;          use Multprec_Integer_Numbers;
with Multprec_Integer_Numbers_io;       use Multprec_Integer_Numbers_io;
with Multprec_Floating_Numbers;         use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;      use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers;          use Multprec_Complex_Numbers;
with Multprec_Complex_Numbers_io;       use Multprec_Complex_Numbers_io;
with Multprec_Random_Numbers;           use Multprec_Random_Numbers;
with Extended_Random_Numbers;           use Extended_Random_Numbers;

procedure ts_extran is

-- DESCRIPTION :
--   Test on creating random extensions of numbers.

  procedure Test_Natural_Random_Extensions is

    size1,size2 : natural32;
    n1,n2 : Natural_Number;
    ans : character;

  begin
    loop
      put("Give size of first number : "); get(size1);
      n1 := Random(size1);
      put("The 1st random : "); put(n1); new_line;
      put("Give the size of second number : "); get(size2);
      n2 := Extended_Random(n1,size2);
      put("The 2nd random : "); put(n2); new_line;
      put("Do you want more tests ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      Clear(n1); Clear(n2);
    end loop;
  end Test_Natural_Random_Extensions;

  procedure Test_Integer_Random_Extensions is

    size1,size2 : natural32;
    n1,n2 : Integer_Number;
    ans : character;

  begin
    loop
      put("Give size of first number : "); get(size1);
      n1 := Random(size1);
      put("The 1st random : "); put(n1); new_line;
      put("Give the size of second number : "); get(size2);
      n2 := Extended_Random(n1,size2);
      put("The 2nd random : "); put(n2); new_line;
      put("Do you want more tests ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      Clear(n1); Clear(n2);
    end loop;
  end Test_Integer_Random_Extensions;

  procedure Test_Floating_Random_Extensions is

    size1,size2 : natural32;
    n1,n2 : Floating_Number;
    ans : character;

  begin
    loop
      put("Give size of first number : "); get(size1);
      n1 := Random(size1);
      put("The 1st random : "); put(n1); new_line;
      put("Give the size of second number : "); get(size2);
      n2 := Extended_Random(n1,size2);
      put("The 2nd random : "); put(n2); new_line;
      put("Do you want more tests ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      Clear(n1); Clear(n2);
    end loop;
  end Test_Floating_Random_Extensions;

  procedure Test_Complex_Random_Extensions is

    size1,size2 : natural32;
    n1,n2 : Complex_Number;
    ans : character;

  begin
    loop
      put("Give size of first number : "); get(size1);
      n1 := Random(size1);
      put("The 1st random : "); put(n1); new_line;
      put("Give the size of second number : "); get(size2);
      n2 := Extended_Random(n1,size2);
      put("The 2nd random : "); put(n2); new_line;
      put("Do you want more tests ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      Clear(n1); Clear(n2);
    end loop;
  end Test_Complex_Random_Extensions;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Generating random extensions of multi-precision numbers.");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line("  0. exit this program;");
      put_line("  1. test random extension of natural numbers;");
      put_line("  2. test random extension of integer numbers;");
      put_line("  3. test random extension of floating numbers.");
      put_line("  4. test random extension of complex numbers.");
      put("Type 0, 1, 2, 3 or 4 to make your choice : ");
      Ask_Alternative(ans,"01234");
      exit when (ans = '0');
      new_line;
      case ans is
        when '1' => Test_Natural_Random_Extensions;
        when '2' => Test_Integer_Random_Extensions;
        when '3' => Test_Floating_Random_Extensions;
        when '4' => Test_Complex_Random_Extensions;
        when others => null;
      end case;
    end loop;
  end Main;

begin
  Main;
end ts_extran;
