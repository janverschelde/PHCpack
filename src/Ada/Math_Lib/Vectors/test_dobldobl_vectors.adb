with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Double_Double_Vectors;
with Double_Double_Vectors_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;
with DoblDobl_Random_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;
with Multprec_Complex_Vectors;
with Multprec_Complex_Vectors_io;
with DoblDobl_Complex_Vectors_cv;        use DoblDobl_Complex_Vectors_cv;

package body Test_DoblDobl_Vectors is

  procedure Test_Double_Double_Vectors is

    use Double_Double_Vectors;
    use Double_Double_Vectors_io;

    n : integer32 := 0;
    m : natural32 := 0;

  begin
    new_line;
    put("Give the dimension : "); get(n);
    put("Give the magnitude : "); get(m);
    declare
      v : Vector(1..n);
    begin
      if m = 1
       then v := DoblDobl_Random_Vectors.Random_Vector(1,n);
       else v := DoblDobl_Random_Vectors.Random_Vector(1,n,m);
      end if;
      put_line("a random vector : "); put_line(v);
    end;
  end Test_Double_Double_Vectors;

  procedure Test_Double_Double_Complex_Vectors is

    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Vectors_io;

    n : integer32 := 0;
    m : natural32 := 0;

  begin
    new_line;
    put("Give the dimension : "); get(n);
    put("Give the magnitude : "); get(m);
    declare
      v : Vector(1..n);
    begin
      if m = 1 
       then v := DoblDobl_Random_Vectors.Random_Vector(1,n);
       else v := DoblDobl_Random_Vectors.Random_Vector(1,n,m);
      end if;
      put_line("a random vector : "); put_line(v);
    end;
  end Test_Double_Double_Complex_Vectors;

  procedure Test_Convertors is

    use Standard_Complex_Vectors_io;
    use DoblDobl_Complex_Vectors_io;
    use Multprec_Complex_Vectors_io;

    n : integer32 := 0;

  begin
    new_line; put("Give the dimension : "); get(n);
    declare
      dd_v : constant DoblDobl_Complex_Vectors.Vector(1..n)
           := DoblDobl_Random_Vectors.Random_Vector(1,n);
      st_v : constant Standard_Complex_Vectors.Vector(dd_v'range)
           := DoblDobl_Complex_to_Standard(dd_v);
      mp_v : constant Multprec_Complex_Vectors.Vector(dd_v'range)
           := DoblDobl_Complex_to_Multprec(dd_v);
    begin
      put_line("a random double double vector : "); put_line(dd_v);
      put_line("-> as standard doubles : "); put_line(st_v);
      put_line("-> as multiprecision : "); put_line(mp_v);
    end;
  end Test_Convertors;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Interactive testing of vectors of double doubles.");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line("  0. exit this program;");
      put_line("  1. random vector of double double numbers;");
      put_line("  2. random vector of complex double double numbers;");
      put_line("  3. convert complex vectors.");
      put("Make your choice (0, 1, 2, or 3) : ");
      Ask_Alternative(ans,"0123");
      exit when (ans = '0');
      case ans is
        when '1' => Test_Double_Double_Vectors;
        when '2' => Test_Double_Double_Complex_Vectors;
        when '3' => Test_Convertors;
        when others => null;
      end case;
    end loop;
  end Main;

end Test_DoblDobl_Vectors;
