with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Deca_Double_Numbers;               use Deca_Double_Numbers;
with DecaDobl_Complex_Numbers;
with DecaDobl_Complex_Vectors_io;
with DecaDobl_Complex_Series;
with DecaDobl_Complex_Series_io;        use DecaDobl_Complex_Series_io;
with DecaDobl_Complex_Random_Series;

package body Test_DecaDobl_Complex_Series is

  procedure DecaDobl_Construct is

    i : constant integer := 123;
    first : constant DecaDobl_Complex_Series.Series
          := DecaDobl_Complex_Series.Create(i);
    second : DecaDobl_Complex_Series.Series(4);
    third : DecaDobl_Complex_Series.Series(4);

  begin
    text_io.put_line("The first series, the constant 123 : ");
    DecaDobl_Complex_Vectors_io.put_line(first.cff);
    text_io.put_line("The second series : ");
    second.cff(0) := DecaDobl_Complex_Numbers.Create(integer32(1));
    second.cff(1) := DecaDobl_Complex_Numbers.Create(integer32(-1));
    second.cff(2) := DecaDobl_Complex_Numbers.Create(integer32(0));
    second.cff(3) := DecaDobl_Complex_Numbers.Create(integer32(0));
    second.cff(4) := DecaDobl_Complex_Numbers.Create(integer32(0));
    DecaDobl_Complex_Vectors_io.put_line(second.cff);
    third := DecaDobl_Complex_Series.Inverse(second);
    text_io.put_line("The inverse of the second series : ");
    DecaDobl_Complex_Vectors_io.put_line(third.cff);
  end DecaDobl_Construct;

  procedure DecaDobl_Test_Creation ( degree : in integer32 ) is

    use DecaDobl_Complex_Numbers;
    use DecaDobl_Complex_Series;

    s : constant Link_to_Series := Create(1,degree);
    t,x,y,z : Link_to_Series;
    minone : constant deca_double := create(-1.0);

  begin
    put("One as series of degree "); put(degree,1); put_line(" :");
    put(s);
    t := Create(1,degree);
    t.cff(1) := Create(minone);
    put_line("The series 1 - t :"); put(t); 
    x := s/t;
    put_line("The series 1/(1-t) : "); put(x);
    y := x*t;
    put_line("Verifying multiplication with inverse : "); put(y);
    z := t*x;
    put_line("Verifying commutativity : "); put(z);
  end DecaDobl_Test_Creation;

  procedure DecaDobl_Test_Arithmetic ( degree : in integer32 ) is

    use DecaDobl_Complex_Series;

    a,b,c : Series(degree);
    ans : character;

  begin
    a := DecaDobl_Complex_Random_Series.Random_Series(degree);
    put_line("The first random series A :"); put(a);
    b := DecaDobl_Complex_Random_Series.Random_Series(degree);
    put_line("The second random series B :"); put(b);
    c := a+b;
    put_line("The sum A + B :"); put(c);
    c := c-a;
    put_line("The sum A + B - A :"); put(c); 
    new_line;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      c := a*b;
      put_line("The product A*B :"); put(c);
      c := c/a;
      put_line("The product A*B/A :"); put(c);
    end if;
  end DecaDobl_Test_Arithmetic;

  procedure Main is

    ans : character;
    degree : integer32 := 0;

  begin
    new_line;
    put_line("MENU with testing operations :");
    put_line("  0. test the basic construct methods");
    put_line("  1. test the computation of 1/(1-t) for any degree");
    put_line("  2. test arithmetic");
    put("Type 0, 1, or 2 to select a test : ");
    Ask_Alternative(ans,"012");
    if ans /= '0' then
      new_line;
      put("Give the degree of the series : "); get(degree);
    end if;
    new_line;
    case ans is 
      when '0' => DecaDobl_Construct;
      when '1' => DecaDobl_Test_Creation(degree);
      when '2' => DecaDobl_Test_Arithmetic(degree);
      when others => null;
    end case;
  end Main;

end Test_DecaDobl_Complex_Series;
