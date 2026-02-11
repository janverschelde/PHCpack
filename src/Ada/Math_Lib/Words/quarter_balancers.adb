with Ada.Text_IO;                        use Ada.Text_IO;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Bits_of_Doubles;

package body Quarter_Balancers is

  function Is_Balanced
             ( x,y : double_float; threshold : integer32 := 13;
               vrblvl : integer32 := 0 ) return boolean is

    res : boolean;
    ex : constant integer32 := double_float'exponent(x);
    ey : constant integer32 := double_float'exponent(y);
    dxy : constant integer32 := ex - ey;

  begin
    if vrblvl > 0 then
      put("-> in Test_Quarter_Balancers.is_balanced, threshold : ");
      put(threshold,1); put_line(" ...");
    end if;
    if vrblvl > 0 then
      put("x : "); put(x); put(" has exponent "); put(ex,1); new_line;
      put("y : "); put(y); put(" has exponent "); put(ey,1); new_line;
      put("difference of exponents : "); put(dxy,1);
    end if;
    res := not (dxy > threshold);
    if vrblvl > 0 then
      if res
       then put_line(" balanced");
       else put_line(" unbalanced");
      end if;
    end if;
    return res;
  end Is_Balanced;

  function Is_Quarter_Balanced
             ( x,y : double_float; vrblvl : integer32 := 0 )
             return boolean is
  begin
    if vrblvl > 0
     then put_line("-> in Test_Quarter_Balancers.is_quarter_balanced ...");
    end if;
    return Is_Balanced(x,y,vrblvl=>vrblvl);
  end Is_Quarter_Balanced;

  procedure Balance ( x,y : in out double_float;
                      threshold : in integer32 := 13;
                      vrblvl : in integer32 := 0 ) is

    exn : constant integer32 := integer32(double_float'exponent(x));
    bit : constant double_float
        := double_float'compose(1.0, exn - threshold);

  begin
    if vrblvl > 0 then
      put("-> in Test_Quarter_Balancers.balance, threshold : ");
      put(threshold,1); put_line(" ...");
    end if;
    if x = 0.0 then
      if vrblvl > 0
       then put_line("x is zero, no point to balance.");
      end if;
      return;
    end if;
    if y = 0.0 then
      if vrblvl > 0
       then put_line("y is zero, no point to balance.");
      end if;
      return;
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
  end Balance;

  procedure Quarter_Balance
              ( x,y : in out double_float; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in Test_Quarter_Balancers.quarter_balance 0 ...");
    end if;
    Balance(x,y,vrblvl=>vrblvl);
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

  procedure Octo_Balance
              ( x0,x1,x2,x3,x4,x5,x6,x7 : in out double_float;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in Test_Quarter_Balancers.octo_balance ...");
    end if;
    if not Is_Balanced(x0,x1,7,vrblvl-1)
     then Balance(x0,x1,7,vrblvl-1);
    end if;
    if not Is_Balanced(x1,x2,7,vrblvl-1)
     then Balance(x1,x2,7,vrblvl-1);
    end if;
    if not Is_Balanced(x2,x3,7,vrblvl-1)
     then Balance(x2,x3,7,vrblvl-1);
    end if;
    if not Is_Balanced(x3,x4,7,vrblvl-1)
     then Balance(x3,x4,7,vrblvl-1);
    end if;
    if not Is_Balanced(x4,x5,6,vrblvl-1)
     then Balance(x4,x5,6,vrblvl-1);
    end if;
    if not Is_Balanced(x5,x6,6,vrblvl-1)
     then Balance(x5,x6,6,vrblvl-1);
    end if;
    if not Is_Balanced(x6,x7,6,vrblvl-1)
     then Balance(x6,x7,6,vrblvl-1);
    end if;
  end Octo_Balance;

end Quarter_Balancers;
