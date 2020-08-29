with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Octo_Double_Numbers_io;             use Octo_Double_Numbers_io;
with OctoDobl_Complex_Numbers;           use OctoDobl_Complex_Numbers;
with OctoDobl_Complex_Numbers_io;        use OctoDobl_Complex_Numbers_io;
with OctoDobl_Random_Numbers;

package body Test_OctoDobl_Complex is

  procedure Test_io is

    c : Complex_Number;

  begin
    put("Give a complex number : "); get(c);
    put_line("-> the real part : "); put(REAL_PART(c)); new_line;
    put_line("-> the imaginary part : "); put(IMAG_PART(c)); new_line;
    put_line("-> your number :"); put(c); new_line;
  end Test_io;

  procedure Test_Addition_and_Subtraction is

    x : constant Complex_Number := OctoDobl_Random_Numbers.Random;
    y : constant Complex_Number := OctoDobl_Random_Numbers.Random;
    s,d : Complex_Number;

  begin
    new_line;
    put_line("Testing x + y - x for random x and y ...");
    put_line("x = "); put(x); new_line;
    put_line("y = "); put(y); new_line;
    s := x + y;
    put_line("x + y = "); put(s); new_line;
    d := s - x;
    put_line("(x + y) - x = "); put(d); new_line;
  end Test_Addition_and_Subtraction;

  procedure Test_Multiplication_and_Division is

    x : constant Complex_Number := OctoDobl_Random_Numbers.Random;
    y : constant Complex_Number := OctoDobl_Random_Numbers.Random;
    p,q : Complex_Number;

  begin
    new_line;
    put_line("Testing x * y / x for random x and y ...");
    put_line("x = "); put(x); new_line;
    put_line("y = "); put(y); new_line;
    p := x * y;
    put_line("x * y = "); put(p); new_line;
    q := p / x;
    put_line("x * y / x = "); put(q); new_line;
  end Test_Multiplication_and_Division;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Testing octo double complex arithmetic ...");
    put_line("  1. test input and output");
    put_line("  2. test addition and subtraction");
    put_line("  3. test multiplication and division");
    put("Type 1, 2, or 3 to select a test : ");
    Ask_Alternative(ans,"123");
    case ans is
      when '1' => Test_io;
      when '2' => Test_Addition_and_Subtraction;
      when '3' => Test_Multiplication_and_Division;
      when others => null;
    end case;
  end Main;

end Test_OctoDobl_Complex;
