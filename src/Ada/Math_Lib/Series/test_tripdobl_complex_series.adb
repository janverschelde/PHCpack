with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Triple_Double_Numbers;             use Triple_Double_Numbers;
with TripDobl_Complex_Numbers;
with TripDobl_Complex_Vectors_io;
with TripDobl_Complex_Series;
with TripDobl_Complex_Series_io;        use TripDobl_Complex_Series_io;
with TripDobl_Complex_Random_Series;

package body Test_TripDobl_Complex_Series is

  procedure TripDobl_Construct is

    i : constant integer := 123;
    first : constant TripDobl_Complex_Series.Series
          := TripDobl_Complex_Series.Create(i);
    second : TripDobl_Complex_Series.Series(4);
    third : TripDobl_Complex_Series.Series(4);

  begin
    text_io.put_line("The first series, the constant 123 : ");
    TripDobl_Complex_Vectors_io.put_line(first.cff);
    text_io.put_line("The second series : ");
    second.cff(0) := TripDobl_Complex_Numbers.Create(integer32(1));
    second.cff(1) := TripDobl_Complex_Numbers.Create(integer32(-1));
    second.cff(2) := TripDobl_Complex_Numbers.Create(integer32(0));
    second.cff(3) := TripDobl_Complex_Numbers.Create(integer32(0));
    second.cff(4) := TripDobl_Complex_Numbers.Create(integer32(0));
    TripDobl_Complex_Vectors_io.put_line(second.cff);
    third := TripDobl_Complex_Series.Inverse(second);
    text_io.put_line("The inverse of the second series : ");
    TripDobl_Complex_Vectors_io.put_line(third.cff);
  end TripDobl_Construct;

  procedure TripDobl_Test_Creation ( degree : in integer32 ) is

    use TripDobl_Complex_Numbers;
    use TripDobl_Complex_Series;

    s : constant Link_to_Series := Create(1,degree);
    t,x,y,z : Link_to_Series;
    minone : constant triple_double := create(-1.0);

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
  end TripDobl_Test_Creation;

  procedure TripDobl_Test_Arithmetic ( degree : in integer32 ) is

    use TripDobl_Complex_Series;

    a,b,c : Series(degree);
    ans : character;

  begin
    a := TripDobl_Complex_Random_Series.Random_Series(degree);
    put_line("The first random series A :"); put(a);
    b := TripDobl_Complex_Random_Series.Random_Series(degree);
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
  end TripDobl_Test_Arithmetic;

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
      when '0' => TripDobl_Construct;
      when '1' => TripDobl_Test_Creation(degree);
      when '2' => TripDobl_Test_Arithmetic(degree);
      when others => null;
    end case;
  end Main;

end Test_TripDobl_Complex_Series;
