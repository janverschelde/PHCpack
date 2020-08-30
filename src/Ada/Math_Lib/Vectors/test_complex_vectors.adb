with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;
with Standard_Random_Vectors;
with Multprec_Complex_Vectors;
with Multprec_Complex_Vectors_io;

package body Test_Complex_Vectors is

  procedure Test_Standard_Vectors_io is

    use Standard_Complex_Vectors;
    use Standard_Complex_Vectors_io;

    n : integer32 := 0;

  begin
    put("Give the dimension : "); get(n);
    declare
      nv : Vector(1..n);
    begin
      put("Give "); put(n,1); put_line(" complex numbers : "); get(nv);
      put_line("Your vector : "); put_line(nv);
    end;
  end Test_Standard_Vectors_io;

  procedure Test_Multprec_Vectors_io is

    use Multprec_Complex_Vectors;
    use Multprec_Complex_Vectors_io;

    n : integer32 := 0;

  begin
    put("Give the dimension : "); get(n);
    declare
      nv : Vector(1..n);
    begin
      put("Give "); put(n,1); put_line(" complex numbers : "); get(nv);
      put_line("Your vector : "); put_line(nv);
    end;
  end Test_Multprec_Vectors_io;

  procedure Test_Random_Vectors is

    n : integer32 := 0;
    m : natural32 := 0;

  begin
    put("Give the dimension : "); get(n);
    put("Give the magnitude : "); get(m);
    declare
      v : constant Standard_Complex_Vectors.Vector(1..n)
        := Standard_Random_Vectors.Random_Vector(1,n,m);
    begin
      put_line("The random vector : ");
      Standard_Complex_Vectors_io.put_line(v);
    end;
  end Test_Random_Vectors;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Interactive testing of vectors of complex numbers.");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line("  0. exit this program.");
      put_line("  1. io of vectors of standard complex numbers;");
      put_line("  2. io of vectors of multi-precision complex numbers;");
      put_line("  3. generate a random standard complex vector.");
      put("Make your choice (0, 1, 2, or 3) : ");
      Ask_Alternative(ans,"0123");
      exit when (ans = '0');
      case ans is
        when '1' => Test_Standard_Vectors_io;
        when '2' => Test_Multprec_Vectors_io;
        when '3' => Test_Random_Vectors;
        when others => null;
      end case;
    end loop;
  end Main;

end Test_Complex_Vectors;
