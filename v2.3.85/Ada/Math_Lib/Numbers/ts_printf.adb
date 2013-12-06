with text_io;                            use text_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;

procedure ts_printf is

-- DESCRIPTION :
--   Test to develop a nice way to print out floating point numbers
--   in the scientific format with little e, as done in C.

-- INTERFACE with C

--  procedure printf ( s : in string; f : in double_float );
--  pragma Import(C,printf,"printf");

--  procedure scanf ( s : in string; f : in out double_float );
--  pragma Import(C,scanf,"scanf");

-- Main program :

  procedure Main is

    f : double_float := 0.0;

  begin

    put("Give a float in Ada format : "); get(f);
    put("float in Ada format : "); put(f); new_line;
    put("float in  C  format : "); printf("%.14e",f); new_line;
    put("Give float in scientific C format : "); 
    scanf("%lf",f);
    put("float in Ada format : "); put(f); new_line;
    put("float in  C  format : "); printf("%.14e",f); new_line;

  end Main;

begin
  Main;
end ts_printf;
