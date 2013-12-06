with text_io;                            use text_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;

procedure ts_cprintf is

-- DESCRIPTION :
--   Test to develop a nice way to print out floating point numbers
--   in the scientific format with little e, as done in C.

  f : double_float;

begin
  put("Give a float : "); get(f);
  put("your float   : "); put(f); new_line;
end ts_cprintf;
