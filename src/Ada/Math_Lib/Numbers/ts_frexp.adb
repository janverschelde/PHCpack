with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;

procedure ts_frexp is

-- DESCRIPTION :
--   Test to get the fraction and the exponent from a double,
--   the equivalent constructions of frexp() and ldexp() in C.

  x : constant double_float := 3.141592653589793;
  f : constant double_float := double_float'fraction(x);
  e : constant integer32 := integer32(double_float'exponent(x));
  c : constant double_float := double_float'compose(f, e);
  s : constant double_float := double_float'compose(f, 52);
  m : constant integer64 := integer64(double_float'truncation(s));

begin
  put("  The 64-bit float : ");
  put("x : "); put(x); new_line;
  put("->    the fraction : ");
  put("f : "); put(f); new_line;
  put("->    the exponent :");
  put("e : "); put(e); new_line;
  put("-> composed number : ");
  put("c : "); put(c); new_line;
  put("-> 52-bit fraction : ");
  put("m : "); put(m); new_line;
end ts_frexp;
