with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Quad_Double_Vectors;
with Quad_Double_Vectors_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;
with QuadDobl_Random_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;
with Multprec_Complex_Vectors;
with Multprec_Complex_Vectors_io;
with QuadDobl_Complex_Vectors_cv;        use QuadDobl_Complex_Vectors_cv;

package body Test_QuadDobl_Vectors is

  procedure Test_Quad_Double_Vectors is

    use Quad_Double_Vectors;
    use Quad_Double_Vectors_io;

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
       then v := QuadDobl_Random_Vectors.Random_Vector(1,n);
       else v := QuadDobl_Random_Vectors.Random_Vector(1,n,m);
      end if;
      put_line("a random vector : "); put_line(v);
    end;
  end Test_Quad_Double_Vectors;

  procedure Test_Quad_Double_Complex_Vectors is

    use QuadDobl_Complex_Vectors;
    use QuadDobl_Complex_Vectors_io;

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
       then v := QuadDobl_Random_Vectors.Random_Vector(1,n);
       else v := QuadDobl_Random_Vectors.Random_Vector(1,n,m);
      end if;
      put_line("a random vector : "); put_line(v);
    end;
  end Test_Quad_Double_Complex_Vectors;

  procedure Test_Convertors is

    use Standard_Complex_Vectors_io;
    use QuadDobl_Complex_Vectors_io;
    use Multprec_Complex_Vectors_io;

    n : integer32 := 0;

  begin
    new_line; put("Give the dimension : "); get(n);
    declare
      qd_v : constant QuadDobl_Complex_Vectors.Vector(1..n)
           := QuadDobl_Random_Vectors.Random_Vector(1,n);
      st_v : constant Standard_Complex_Vectors.Vector(qd_v'range)
           := QuadDobl_Complex_to_Standard(qd_v);
      mp_v : constant Multprec_Complex_Vectors.Vector(qd_v'range)
           := QuadDobl_Complex_to_Multprec(qd_v);
    begin
      put_line("a random quad double vector : "); put_line(qd_v);
      put_line("-> as standard doubles : "); put_line(st_v);
      put_line("-> as multiprecision : "); put_line(mp_v);
    end;
  end Test_Convertors;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Interactive testing of vectors of quad doubles.");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line("  0. exit this program;");
      put_line("  1. random vector of quad double numbers;");
      put_line("  2. random vector of complex quad double numbers;");
      put_line("  3. convert complex vectors.");
      put("Make your choice (0, 1, 2, or 3) : ");
      Ask_Alternative(ans,"0123");
      exit when (ans = '0');
      case ans is
        when '1' => Test_Quad_Double_Vectors;
        when '2' => Test_Quad_Double_Complex_Vectors;
        when '3' => Test_Convertors;
        when others => null;
      end case;
    end loop;
  end Main;

end Test_QuadDobl_Vectors;
