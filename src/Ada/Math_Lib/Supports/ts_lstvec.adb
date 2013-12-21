with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_Vectors;           use Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Lists_of_Integer_Vectors_io;        use Lists_of_Integer_Vectors_io;
with Lists_of_Floating_Vectors;          use Lists_of_Floating_Vectors;
with Lists_of_Floating_Vectors_io;       use Lists_of_Floating_Vectors_io;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Arrays_of_Integer_Vector_Lists_io;  use Arrays_of_Integer_Vector_Lists_io;
with Arrays_of_Floating_Vector_Lists;    use Arrays_of_Floating_Vector_Lists;
with Arrays_of_Floating_Vector_Lists_io; use Arrays_of_Floating_Vector_Lists_io;

procedure ts_lstvec is

-- DESCRIPTION :
--   Test on operations on lists of links to vectors.

  procedure Test_Integer_List_io is

    m,n : natural32 := 0;
    l : Lists_of_Integer_Vectors.List;

  begin
    new_line;
    put_line("Testing input/output of lists of links to integer vectors.");
    new_line;
    put("Give the dimension of the vectors : "); get(n);
    put("Give the number of vectors : "); get(m);
    put("Give "); put(m,1); put(" "); put(n,1); put_line("-vectors :");
    get(n,m,l);
    put_line("Your list : "); put(l);
  end Test_Integer_List_io;

  procedure Test_Floating_List_io is

    m,n : natural32 := 0;
    l : Lists_of_Floating_Vectors.List;

  begin
    new_line;
    put_line("Testing input/output of lists of links to floating vectors.");
    new_line;
    put("Give the dimension of the vectors : "); get(n);
    put("Give the number of vectors : "); get(m);
    put("Give "); put(m,1); put(" "); put(n,1); put_line("-vectors :");
    get(n,m,l);
    put_line("Your list : "); put(l);
  end Test_Floating_List_io;

  procedure Test_Integer_Array_List_io is

    n : natural32 := 0;

  begin
    new_line;
    put_line("Testing input/output of arrays of lists of integer vectors.");
    new_line;
    put("Give the dimension : "); get(n);
    declare
      l : Arrays_of_Integer_Vector_Lists.Array_of_Lists(1..integer32(n));
      m : Standard_Natural_Vectors.Vector(1..integer32(n));
    begin
      put("Give the cardinalities of the lists : "); get(m);
      put_line("Give the lists : "); get(n,m,l);
      put_line("Your lists : "); put(l);
    end;
  end Test_Integer_Array_List_io;

  procedure Test_Floating_Array_List_io is

    n : natural32 := 0;

  begin
    new_line;
    put_line("Testing input/output of arrays of lists of floating vectors.");
    new_line;
    put("Give the dimension : "); get(n);
    declare
      l : Arrays_of_Floating_Vector_Lists.Array_of_Lists(1..integer32(n));
      m : Standard_Natural_Vectors.Vector(1..integer32(n));
    begin
      put("Give the cardinalities of the lists : "); get(m);
      put_line("Give the lists : "); get(n,m,l);
      put_line("Your lists : "); put(l);
    end;
  end Test_Floating_Array_List_io;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Interactive testing of lists of links to vectors.");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line("  0. Exit this program.");
      put_line("  1. i/o for lists of links to standard integer vectors.");
      put_line("  2. i/o for lists of links to standard floating vectors.");
      put_line("  3. i/o for arrays of standard integer vector lists.");
      put_line("  4. i/o for arrays of standard floating vector lists.");
      put("Type 0,1,2,3 or 4 to select : "); get(ans);
      exit when (ans = '0');
      case ans is
        when '1' => Test_Integer_List_io;
        when '2' => Test_Floating_List_io;
        when '3' => Test_Integer_Array_List_io;
        when '4' => Test_Floating_Array_List_io;
        when others => null;
      end case;
    end loop;
  end Main;

begin
  Main;
end ts_lstvec;
