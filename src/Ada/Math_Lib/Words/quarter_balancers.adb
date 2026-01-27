with Ada.Text_IO;                        use Ada.Text_IO;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Bits_of_Doubles;

package body Quarter_Balancers is

  function Is_Quarter_Balanced
             ( x,y : double_float; vrblvl : integer32 := 0 )
             return boolean is

    res : boolean;
    ex : constant integer32 := double_float'exponent(x);
    ey : constant integer32 := double_float'exponent(y);
    dxy : constant integer32 := ex - ey;

  begin
    if vrblvl > 0
     then put_line("-> in Test_Quarter_Balancers.is_quarter_balanced ...");
    end if;
    if vrblvl > 0 then
      put("x : "); put(x); put(" has exponent "); put(ex,1); new_line;
      put("y : "); put(y); put(" has exponent "); put(ey,1); new_line;
      put("difference of exponents : "); put(dxy,1);
    end if;
    res := not (dxy > 13);
    if res
     then put_line(" balanced");
     else put_line(" unbalanced");
    end if;
    return res;
  end Is_Quarter_Balanced;

  procedure Quarter_Balance
              ( x,y : in out double_float; vrblvl : in integer32 := 0 ) is

    exn : constant integer32 := integer32(double_float'exponent(x));
    bit : constant double_float := double_float'compose(1.0, exn - 13);

  begin
    if vrblvl > 0
     then put_line("-> in Test_Quarter_Balancers.quarter_balance 0 ...");
    end if;
    if vrblvl > 0 then
      put("b x : "); Bits_of_Doubles.write_52bits_expo(x); new_line;
      put("bit : "); Bits_of_Doubles.write_52bits_expo(bit); new_line;
      put("b y : "); Bits_of_Doubles.write_52bits_expo(y); new_line;
    end if;
    x := x - bit;
    y := y + bit;
    if vrblvl > 0 then
      put("b x : "); Bits_of_Doubles.write_52bits_expo(x); new_line;
      put("b y : "); Bits_of_Doubles.write_52bits_expo(y); new_line;
    end if;
  end Quarter_Balance;

  procedure Quarter_Balance
              ( x0,x1,x2,x3 : in out double_float;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in Test_Quarter_Balancers.quarter_balance 1 ...");
    end if;
    if not Is_Quarter_Balanced(x0,x1,vrblvl-1)
     then Quarter_Balance(x0,x1,vrblvl-1);
    end if;
    if not Is_Quarter_Balanced(x1,x2,vrblvl-1)
     then Quarter_Balance(x1,x2,vrblvl-1);
    end if;
    if not Is_Quarter_Balanced(x2,x3,vrblvl-1)
     then Quarter_Balance(x2,x3,vrblvl-1);
    end if;
  end Quarter_Balance;

end Quarter_Balancers;
