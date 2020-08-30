with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;
with Standard_Integer64_Vectors;
with Standard_Integer64_Vectors_io;
with Standard_Integer_VecVecs;
with Standard_Integer_VecVecs_io;
with Standard_Integer64_VecVecs;
with Standard_Integer64_VecVecs_io;
with Multprec_Integer_Vectors;
with Multprec_Integer_Vectors_io;
with Multprec_Integer_Vecvecs;
with Multprec_Integer_Vecvecs_io;

package body Test_Integer_Vectors is

  procedure Test_Standard_Vectors_io is

    use Standard_Integer_Vectors;
    use Standard_Integer_Vectors_io;

    n : integer32 := 0;

  begin
    put("Give the dimension : "); get(n);
    declare
      nv : Vector(1..n);
    begin
      put("Give "); put(n,1); put_line(" integer numbers : "); get(nv);
      put("Your vector : "); put(nv); new_line;
    end;
  end Test_Standard_Vectors_io;

  procedure Test_Standard64_Vectors_io is

    use Standard_Integer64_Vectors;
    use Standard_Integer64_Vectors_io;

    n : integer32 := 0;

  begin
    put("Give the dimension : "); get(n);
    declare
      nv : Vector(1..n);
    begin
      put("Give "); put(n,1); put_line(" integer numbers : "); get(nv);
      put("Your vector : "); put(nv); new_line;
    end;
  end Test_Standard64_Vectors_io;

  procedure Test_Standard_VecVecs_io is

    use Standard_Integer_VecVecs;
    use Standard_Integer_VecVecs_io;

    n : integer32 := 0;

  begin
    put("Give the dimension : "); get(n);
    declare
      nv : VecVec(1..n);
    begin
      put("Give "); put(n,1); put_line(" integer vectors : ");
      get(natural32(n),nv);
      put_line("Your vector : "); put(nv); new_line;
    end;
  end Test_Standard_VecVecs_io;

  procedure Test_Standard64_VecVecs_io is

    use Standard_Integer64_VecVecs;
    use Standard_Integer64_VecVecs_io;

    n : integer32 := 0;

  begin
    put("Give the dimension : "); get(n);
    declare
      nv : VecVec(1..n);
    begin
      put("Give "); put(n,1); put_line(" integer vectors : ");
      get(natural32(n),nv);
      put_line("Your vector : "); put(nv); new_line;
    end;
  end Test_Standard64_VecVecs_io;

  procedure Test_Multprec_Vectors_io is

    use Multprec_Integer_Vectors,Multprec_Integer_Vectors_io;

    n : integer32 := 0;

  begin
    put("Give the dimension : "); get(n);
    declare
      nv : Vector(1..n);
    begin
      put("Give "); put(n,1); put_line(" integer numbers : "); get(nv);
      put_line("Your vector : "); put_line(nv);
    end;
  end Test_Multprec_Vectors_io;

  procedure Test_Multprec_VecVecs_io is

    use Multprec_Integer_VecVecs;
    use Multprec_Integer_VecVecs_io;

    n : integer32 := 0;

  begin
    put("Give the dimension : "); get(n);
    declare
      nv : VecVec(1..n);
    begin
      put("Give "); put(n,1); put_line(" integer vectors : ");
      get(natural32(n),nv);
      put_line("Your vector : "); put_line(nv);
    end;
  end Test_Multprec_VecVecs_io;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Interactive testing of vectors of integer numbers.");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line("  0. exit this program.");
      put_line("  1. io of vectors of standard integer numbers.");
      put_line("  2. io of vectors of standard 64-bit integer numbers.");
      put_line("  3. io of vectors of vectors of 64-bit integer numbers.");
      put_line("  4. io of vectors of multi-precision integer numbers.");
      put_line("  5. io of vectors of vectors of multi-precision integers.");
      put("Make your choice (0,1,2,3,4, or 5) : ");
      Ask_Alternative(ans,"0123456");
      new_line;
      case ans is
        when '1' => Test_Standard_Vectors_io;
        when '2' => Test_Standard64_Vectors_io;
        when '3' => Test_Standard_VecVecs_io;
        when '4' => Test_Standard64_VecVecs_io;
        when '5' => Test_Multprec_Vectors_io;
        when '6' => Test_Multprec_VecVecs_io;
        when others => null;
      end case;
      exit when ans = '0';
    end loop;
  end Main;

end Test_Integer_Vectors;
