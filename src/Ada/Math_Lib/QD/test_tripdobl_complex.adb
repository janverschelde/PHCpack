with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Triple_Double_Numbers;              use Triple_Double_Numbers;
with Triple_Double_Numbers_io;           use Triple_Double_Numbers_io;
with TripDobl_Complex_Numbers;           use TripDobl_Complex_Numbers;
with TripDobl_Complex_Numbers_io;        use TripDobl_Complex_Numbers_io;
with TripDobl_Mathematical_Functions;
with TripDobl_Random_Numbers;
with TripDobl_Complex_Numbers_Polar;     use TripDobl_Complex_Numbers_Polar;

package body Test_TripDobl_Complex is

  procedure Test_io is

    c : Complex_Number;

  begin
    put("Give a complex number : "); get(c);
    put_line("-> the real part : "); put(REAL_PART(c)); new_line;
    put_line("-> the imaginary part : "); put(IMAG_PART(c)); new_line;
    put_line("-> your number :"); put(c); new_line;
  end Test_io;

  procedure Test_Addition_and_Subtraction is

    x : constant Complex_Number := TripDobl_Random_Numbers.Random;
    y : constant Complex_Number := TripDobl_Random_Numbers.Random;
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

    x : constant Complex_Number := TripDobl_Random_Numbers.Random;
    y : constant Complex_Number := TripDobl_Random_Numbers.Random;
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

  procedure Test_Random is

    rnc : constant Complex_Number := TripDobl_Random_Numbers.Random1;
    x : constant triple_double := REAL_PART(rnc);
    y : constant triple_double := IMAG_PART(rnc);
    rad : constant triple_double
        := TripDobl_Mathematical_Functions.Radius(x,y);

  begin
    new_line;
    put_line("A random complex number :"); put(rnc); new_line;
    put("Its radius : "); put(rad); new_line;
  end Test_Random;

  procedure Prompt_Complex_Number ( c : out Complex_Number ) is

    re,im : triple_double;

  begin
    put_line("reading a complex number ...");
    skip_line;
    put("Give re : "); get(re);
    put("Give im : "); get(im);
    c := create(re,im);
    put_line("your complex number : "); put(c); new_line;
  end Prompt_Complex_Number;

  procedure Test_Roots is

    d,k : natural32 := 0;
    a,c,prod,err : Complex_Number;
    ans : character;

  begin
    new_line;
    put_line("Solving x^d - c = 0, for a triple double complex number c.");
    new_line;
    put("Give the degree d : "); get(d);
    Prompt_Complex_Number(c);
    loop
      put("Which root do you want ? "); get(k);
      a := Root(c,d,k);
      put("The root is "); put(a); new_line;
      prod := a;
      for j in 2..d loop
        prod := prod*a;
      end loop;
      put("root^d  =   "); put(prod); new_line;
      err := prod - c;
      put("root^d - c = "); put(err,3); new_line;
      put("Do you want other roots ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_Roots;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Testing triple double complex arithmetic ...");
    put_line("  1. test input and output");
    put_line("  2. test addition and subtraction");
    put_line("  3. test multiplication and division");
    put_line("  4. generate a random complex number");
    put_line("  5. test computation of primitive roots");
    put("Type 1, 2, 3, 4, or 5 to select a test : ");
    Ask_Alternative(ans,"12345");
    case ans is
      when '1' => Test_io;
      when '2' => Test_Addition_and_Subtraction;
      when '3' => Test_Multiplication_and_Division;
      when '4' => Test_Random;
      when '5' => Test_Roots;
      when others => null;
    end case;
  end Main;

end Test_TripDobl_Complex;
