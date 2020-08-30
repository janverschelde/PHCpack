with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;
with Standard_Random_Vectors;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecs_io;
with Multprec_Floating_Vectors;
with Multprec_Floating_Vectors_io;
with Multprec_Floating_Vecvecs;
with Multprec_Floating_Vecvecs_io;

package body Test_Floating_Vectors is

  procedure Test_Standard_Vectors_io is

    use Standard_Floating_Vectors;
    use Standard_Floating_Vectors_io;

    n : integer32 := 0;

  begin
    put("Give the dimension : "); get(n);
    declare
      nv : Vector(1..n);
    begin
      put("Give "); put(n,1); put_line(" floating numbers : "); get(nv);
      put("Your vector : "); put(nv); new_line;
    end;
  end Test_Standard_Vectors_io;

  procedure Test_Standard_VecVecs_io is

    use Standard_Floating_VecVecs;
    use Standard_Floating_VecVecs_io;

    n : integer32 := 0;

  begin
    put("Give the dimension : "); get(n);
    declare
      nv : VecVec(1..n);
    begin
      put("Give "); put(n,1); put_line(" floating vectors : ");
      get(natural32(n),nv);
      put_line("Your vector : "); put(nv); new_line;
    end;
  end Test_Standard_VecVecs_io;

  procedure Test_Multprec_Vectors_io is

    use Multprec_Floating_Vectors;
    use Multprec_Floating_Vectors_io;

    n : integer32 := 0;

  begin
    put("Give the dimension : "); get(n);
    declare
      nv : Vector(1..n);
    begin
      put("Give "); put(n,1); put_line(" floating numbers : "); get(nv);
      put_line("Your vector : "); put_line(nv);
    end;
  end Test_Multprec_Vectors_io;

  procedure Test_Multprec_VecVecs_io is

    use Multprec_Floating_VecVecs;
    use Multprec_Floating_VecVecs_io;

    n : integer32 := 0;

  begin
    put("Give the dimension : "); get(n);
    declare
      nv : VecVec(1..n);
    begin
      put("Give "); put(n,1); put_line(" floating vectors : ");
      get(natural32(n),nv);
      put_line("Your vector : "); put_line(nv);
    end;
  end Test_Multprec_VecVecs_io;

  procedure Test_Random_Vectors is

    n : integer32 := 0;
    m : natural32 := 0;

  begin
    put("Give the dimension : "); get(n);
    put("Give the magnitude : "); get(m);
    declare
      v : constant Standard_Floating_Vectors.Vector(1..n)
        := Standard_Random_Vectors.Random_Vector(1,n,m);
    begin
      put_line("The random vector :");
      Standard_Floating_Vectors_io.put_line(v);
    end; 
  end Test_Random_Vectors;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Interactive testing of vectors of floating numbers.");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line("  0. exit this program.");
      put_line("  1. io of vectors of standard floating numbers.");
      put_line("  2. io of vectors of vectors of standard floating numbers.");
      put_line("  3. io of vectors of multi-precision floating numbers.");
      put_line("  4. io of vectors of vectors of multi-precision floats;");
      put_line("  5. generate a random standard floating vector.");
      put("Make your choice (0,1,2,3, 4, or 5) : ");
      Ask_Alternative(ans,"012345");
      case ans is
        when '1' => Test_Standard_Vectors_io;
        when '2' => Test_Standard_VecVecs_io;
        when '3' => Test_Multprec_Vectors_io;
        when '4' => Test_Multprec_VecVecs_io;
        when '5' => Test_Random_Vectors;
        when others => null;
      end case;
      exit when ans = '0';
    end loop;
  end Main;

end Test_Floating_Vectors;
