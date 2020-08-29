with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Communications_with_User;           use Communications_with_User;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with DoblDobl_Random_Numbers;
with DoblDobl_Complex_Numbers_Polar;     use DoblDobl_Complex_Numbers_Polar;

package body Test_DoblDobl_Complex is

--  procedure Write ( c : in complex_number ) is
-- 
--  -- DESCRIPTION :
--  --   Very basic output of a double double complex number c.
--
--  begin
--    write(REAL_PART(c),32); put("  ");
--    write(IMAG_PART(c),32); put("*i"); new_line;
--  end Write;

  procedure Basic_Test is

   -- i : integer;
   -- d : double_double;
    c : complex_number;

  begin
    new_line;
   -- put("Give an integer : "); get(i);
   -- d := Create(i);
   -- c := Create(d,d);
   -- put(c); -- Write(c);
   -- new_line;
   -- skip_line;
    put("Give a complex number c : "); get(c);
    put("-> c : "); put(c); new_line;
  end Basic_Test;

  procedure Test_Addition_and_Subtraction is

    x : constant complex_number := DoblDobl_Random_Numbers.Random;
    y : constant complex_number := DoblDobl_Random_Numbers.Random;
    s,d : complex_number;

  begin
    new_line;
    put_line("Testing x + y - x for random x and y ...");
    put_line("x = "); put(x); new_line;
    put_line("y = "); put(y); new_line;
    s := x + y;
    put_line("x + y = "); put(s); new_line;
    d := s - x;
    put_line("x + y - x = "); put(d); new_line;
  end Test_Addition_and_Subtraction;

  procedure Test_Multiplication_and_Division is

    x : constant complex_number := DoblDobl_Random_Numbers.Random;
    y : constant complex_number := DoblDobl_Random_Numbers.Random;
    p,q : complex_number;

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

  procedure Prompt_Complex_Number ( c : out Complex_Number ) is

    re,im : double_double;

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
    put_line("Solving x^d - c = 0, for a double double complex number c.");
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

  procedure Test_Random is
 
    m : natural32 := 0;
    r : Complex_Number;

  begin
    put("Give the magnitude : "); get(m);
    if m = 1
     then r := DoblDobl_Random_Numbers.Random1;
     else r := DoblDobl_Random_Numbers.Random_Magnitude(m);
    end if;
    put("The random number : "); put(r); new_line;
  end Test_Random;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU for testing double double complex arithmetic :");
    put_line("  0. run a basic test;");
    put_line("  1. test addition and subtraction;");
    put_line("  2. test multiplication and division;");
    put_line("  3. test finding of primitive roots;");
    put_line("  4. generate a random complex number.");
    put("Type 0, 1, 2, 3, or 4 to choose : ");
    Ask_Alternative(ans,"01234");
    case ans is
      when '0' => Basic_Test;
      when '1' => Test_Addition_and_Subtraction;
      when '2' => Test_Multiplication_and_Division;
      when '3' => Test_Roots;
      when '4' => Test_Random;
      when others => null;
    end case;
  end Main;

end Test_DoblDobl_Complex;
