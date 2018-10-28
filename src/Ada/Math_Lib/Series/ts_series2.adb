with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Dense_Series2;
with Standard_Dense_Series2_io;          use Standard_Dense_Series2_io;

procedure ts_series2 is

-- DESCRIPTION :
--   Tests the operations on truncated power series.

  procedure Standard_Test_Creation ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Verifies that 1/(1-t) = 1 + t + t^2 + ...
  --   for a truncated power series with coefficients
  --   in standard double precision.

    use Standard_Complex_Numbers;
    use Standard_Dense_Series2;

    s : constant Series := Create(1,degree);
    t : Series := s;
    x,y,z : Series(degree);

  begin
    put("One as series of degree "); put(degree,1); put_line(" :");
    put(s);
    t.cff(1) := Create(-1.0);
    put_line("The series 1 - t :"); put(t); 
    x := s/t;
    put_line("The series 1/(1-t) : "); put(x);
    y := x*t;
    put_line("Verifying multiplication with inverse : "); put(y);
    z := t*x;
    put_line("Verifying commutativity : "); put(z);
  end Standard_Test_Creation;

  procedure Standard_Test_Division ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the division on random series of the given degree,
  --   in standard double precision.

    use Standard_Dense_Series2;

    a,b,c : Series(degree);

  begin
    put("Give "); put(degree+1,1);
    put_line(" complex numbers for the first series : "); 
    for i in 0..a.deg loop
      get(a.cff(i));
    end loop;
    put_line("The first series : "); put(a);
    new_line;
    put("Give "); put(degree+1,1);
    put_line(" complex numbers for the second series : "); 
    for i in 0..b.deg loop
      get(b.cff(i));
    end loop;
    put_line("The first series : "); put(b);
    c := a/b;
    new_line;
    put_line("The result of the division "); put(c);
  end Standard_Test_Division;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the degree of the series.

    degree : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("MENU with testing operations :");
    put_line("  0. test the computation of 1/(1-t)");
    put_line("  1. test division operation");
    put("Type 0, or 1 to make your choice : ");
    Ask_Alternative(ans,"01");
    new_line;
    put("Give the degree of the series : "); get(degree);
    new_line;
    case ans is
      when '0' => Standard_Test_Creation(degree);
      when '1' => Standard_Test_Division(degree);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_series2;
