with system;
with text_io,integer_io;

procedure ts_longflt is

-- DESCRIPTION :
--   Test on long floats...

  type long_float is digits system.Max_Digits;
  package long_float_io is new text_io.float_io(long_float);
  package long_long_integer_io is new text_io.integer_io(long_long_integer);

 -- type long_long_natural is range 0..long_long_integer'last;

  procedure Main is

    f : long_float;
    use long_float_io;
    use long_long_integer_io;
    use text_io;
    a,b : long_long_integer;

  begin
    put("long_natural'last : ");
    integer_io.put(natural'last); new_line;
    put("integer'last : "); integer_io.put(integer'last); new_line;
    put("long_long_integer'last : ");
    long_long_integer_io.put(long_long_integer'last); new_line;
    put("System.Max_Int : "); put(System.Max_Int,1); new_line;
    put("System.Min_Int : "); put(System.Min_Int,1); new_line;
    put("Give first long long integer : "); get(a);
    put("Give second long long integer : "); get(b);
    put(a,1); put("*"); put(b,1); put(" = "); put(a*b,1); new_line;
    put("System.Max_Digits : "); put(System.Max_Digits,1); new_line;
    put("Give a long float : "); get(f);
    put("  your long float : "); put(f); new_line;
  end Main;

begin
  Main;
end ts_longflt;
