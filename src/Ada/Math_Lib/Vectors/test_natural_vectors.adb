with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;
with Standard_Natural_VecVecs;
with Standard_Natural_VecVecs_io;
with Multprec_Natural_Vectors;
with Multprec_Natural_Vectors_io;
with Multprec_Natural_VecVecs;
with Multprec_Natural_VecVecs_io;

package body Test_Natural_Vectors is

  procedure Test_Standard_Vectors_io is

    use Standard_Natural_Vectors;
    use Standard_Natural_Vectors_io;

    n : integer32 := 0;

  begin
    put("Give the dimension : "); get(n);
    declare
      nv : Vector(1..n);
    begin
      put("Give "); put(n,1); put_line(" natural numbers : "); get(nv);
      put("Your vector : "); put(nv); new_line;
    end;
  end Test_Standard_Vectors_io;

  procedure Test_Standard_Addition is

    use Standard_Natural_Vectors;
    use Standard_Natural_Vectors_io;

    n : integer32 := 0;

  begin
    put("Give the dimension : "); get(n);
    declare
      a,b,sum : Vector(1..n);
    begin
      put("Give "); put(n,1); put_line(" natural numbers : "); get(a);
      put("-> vector a  : "); put(a); new_line;
      put("Give "); put(n,1); put_line(" natural numbers : "); get(b);
      put("-> vector b  : "); put(b); new_line;
      sum := a+b;
      put(" The sum a+b : "); put(sum); new_line;
    end;
  end Test_Standard_Addition;

  procedure Test_Standard_VecVecs_io is

    use Standard_Natural_VecVecs;
    use Standard_Natural_VecVecs_io;

    n : integer32 := 0;

  begin
    put("Give the dimension : "); get(n);
    declare
      nv : VecVec(1..n);
    begin
      put("Give "); put(n,1); put_line(" natural vectors : ");
      get(natural32(n),nv);
      put_line("Your vector : "); put(nv); new_line;
    end;
  end Test_Standard_VecVecs_io;

  procedure Test_Multprec_Vectors_io is

    use Multprec_Natural_Vectors;
    use Multprec_Natural_Vectors_io;

    n : integer32 := 0;

  begin
    put("Give the dimension : "); get(n);
    declare
      nv : Vector(1..n);
    begin
      put("Give "); put(n,1); put_line(" natural numbers : "); get(nv);
      put_line("Your vector : "); put_line(nv);
    end;
  end Test_Multprec_Vectors_io;

  procedure Test_Multprec_VecVecs_io is

    use Multprec_Natural_VecVecs;
    use Multprec_Natural_VecVecs_io;

    n : integer32 := 0;

  begin
    put("Give the dimension : "); get(n);
    declare
      nv : VecVec(1..n);
    begin
      put("Give "); put(n,1); put_line(" natural vectors : ");
      get(natural32(n),nv);
      put_line("Your vector : "); put_line(nv);
    end;
  end Test_Multprec_VecVecs_io;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Interactive testing of vectors of natural numbers.");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line("  0. exit the program.");
      put_line("  1. io of vectors of standard natural numbers.");
      put_line("  2. addition of vectors of standard natural numbers.");
      put_line("  3. io of vectors of vectors of standard natural numbers.");
      put_line("  4. io of vectors of multi-precision natural numbers.");
      put_line("  5. io of vectors of vectors of multi-precision naturals.");
      put("Make your choice (0,1,2,3, or 4) : "); get(ans);
      case ans is
        when '1' => Test_Standard_Vectors_io;
        when '2' => Test_Standard_Addition;
        when '3' => Test_Standard_VecVecs_io;
        when '4' => Test_Multprec_Vectors_io;
        when '5' => Test_Multprec_VecVecs_io;
        when others => null;
      end case;
      exit when (ans = '0');
    end loop;
  end Main;

end Test_Natural_Vectors;
