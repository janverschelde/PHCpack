with text_io;                            use text_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Random_Numbers;
with Bits_of_Doubles;                    use Bits_of_Doubles;

procedure ts_splitdbl is

-- DESCRIPTION :
--   Tests the splitting of a 64-bit floating-point number
--   in four equal parts.

  procedure Main is

  -- DESCRIPTION :
  --   Generates a random 64-bit floating-point number
  --   and splits its fraction in four equal parts.

    x : constant double_float := Standard_Random_Numbers.Random_Magnitude(50);
    x0,x1,x2,x3,s1,s2,s3,err : double_float;
  
  begin
    Split(x,x0,x1,x2,x3);
    s1 := x1 + x0;
    s2 := x2 + s1;
    s3 := x3 + s2;
    err := abs(x-x0);
    put("          x : "); put(x); new_line;
    put("         x0 : "); put(x0);
    put("  error : "); put(err,2); new_line;
    err := abs(x-s1);
    put("          x : "); put(x); new_line;
    put("      x0+x1 : "); put(s1);
    put("  error : "); put(err,2); new_line;
    err := abs(x-s2);
    put("          x : "); put(x); new_line;
    put("   x0+x1+x2 : "); put(s2);
    put("  error : "); put(err,2); new_line;
    err := abs(x-s3);
    put("          x : "); put(x); new_line;
    put("x1+x1+x2+x3 : "); put(s3);
    put("  error : "); put(err,2); new_line;
    put(" b0 : "); write_52bits_expo(x0); new_line;
    put(" b1 : "); write_52bits_expo(x1); new_line;
    put(" b2 : "); write_52bits_expo(x2); new_line;
    put(" b3 : "); write_52bits_expo(x3); new_line;
    put("x : "); write_52bits_expo(x); new_line;
    put("b : "); write_52bits_expo(s3); new_line;
    if Bit_Equal(x,s3)
     then put_line("The sum of the four parts and x are bit equal, okay.");
     else put_line("The sum of the four parts and x are NOT bit equal, bug!");
    end if;
  end Main;

begin
  Main;
end ts_splitdbl;
