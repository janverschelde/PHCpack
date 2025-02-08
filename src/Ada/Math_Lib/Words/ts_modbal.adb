with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Floating_Vectors;
with Standard_Random_Vectors;
with Bits_of_Doubles;                    use Bits_of_Doubles;
with Balanced_Quarter_Doubles;

procedure ts_modbal is

-- DESCRIPTION :
--   Tests the balancing of a double by the addition
--   of a modifier constant.

  function Balancing_Bits return integer64 is

  -- DESCRIPTION :
  --   Returns the sequence 
  --   1 0000 0000 0000 1 0000 0000 0000 1 0000 0000 0000 1 0000 0000 0000
  --   as a 64-bit integer number.

    res : integer64 := 1; -- accumulates the bits

  begin
    for k in 1..4 loop
      for i in 1..12 loop
        res := 2*res;
      end loop;
      exit when (k = 4);
      res := 2*res + 1;
    end loop;
    return res;
  end Balancing_Bits;

  function Modifier_Constant
             ( verbose : in boolean := false ) return double_float is

  -- DESCRIPTION :
  --   Returns the modifier constant to balance a double.

    frc : constant integer64 := Balancing_Bits;
    mrs : constant double_float := double_float(frc);
    res : constant double_float := double_float'compose(mrs,0);

  begin
    if verbose
     then put("b : "); write_52bits_expo(res); new_line;
    end if;
    return res;
  end Modifier_Constant;

  procedure Test_Inner_Product ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates two random vectors of dimension dim
  --   and tests their inner product.

    mct : constant double_float := Modifier_Constant;
    x : constant Standard_Floating_Vectors.Vector(1..dim)
      := Standard_Random_Vectors.Random_Vector(1,dim);
    y : constant Standard_Floating_Vectors.Vector(1..dim)
      := Standard_Random_Vectors.Random_Vector(1,dim);
    xplus,xmin,yplus,ymin : Standard_Floating_Vectors.Vector(1..dim);
    z,zplus,zmin,xsum,ysum : double_float := 0.0;
    result : double_float;

  begin
    for i in 1..dim loop
      z := z + x(i)*y(i);
      xsum := xsum + x(i);
      ysum := ysum + y(i);
    end loop;
    put("z : "); put(z); new_line;
    for i in 1..dim loop
      xplus(i) := x(i) + mct; xmin(i) := x(i) - mct;
      yplus(i) := y(i) + mct; ymin(i) := y(i) - mct;
    end loop;
    for i in 1..dim loop
      zplus := zplus + xplus(i)*yplus(i);
      zmin := zmin + xmin(i)*ymin(i);
    end loop;
    result := (zplus + zmin)/2.0 - double_float(dim)*mct*mct;
    put("z : "); put(result); new_line;
    result := zplus  - mct*xsum - mct*ysum - double_float(dim)*mct*mct;
    put("z : "); put(result); new_line;
  end Test_Inner_Product;

  procedure Test_Balancing ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Checks if a vector of doubles, modified with the constant
  --   is balanced when split.

    mct : constant double_float := Modifier_Constant;
    m0,m1,m2,m3 : double_float;
    x : Standard_Floating_Vectors.Vector(1..dim)
      := Standard_Random_Vectors.Random_Vector(1,dim);
    x0,x1,x2,x3 : double_float;
    xplus,xmin : Standard_Floating_Vectors.Vector(1..dim);
    xbal,xpbal,xmbal,mbal : boolean;
    xcnt,xpcnt,xmcnt : integer32 := 0;

  begin
    Bits_of_Doubles.Split(mct,m0,m1,m2,m3);
    mbal := Balanced_Quarter_Doubles.Is_Balanced(0,m0,m1,m2,m3);
    if mbal
     then put(mct); put_line(" is balanced");
     else put(mct); put_line(" is NOT balanced!");
    end if;
    for i in 1..dim loop
      if x(i) < 0.0
       then x(i) := -x(i);
      end if;
      xplus(i) := x(i) + mct;
      xmin(i) := x(i) - mct;
    end loop;
    for i in 1..dim loop
       Bits_of_Doubles.Split(x(i),x0,x1,x2,x3);
       xbal := Balanced_Quarter_Doubles.Is_Balanced(0,x0,x1,x2,x3);
       Bits_of_Doubles.Split(xplus(i),x0,x1,x2,x3);
       x0 := m0 + x0;
       x1 := m1 + x1;
       x2 := m2 + x2;
       x3 := m3 + x3;
       xpbal := Balanced_Quarter_Doubles.Is_Balanced(0,x0,x1,x2,x3);
       Bits_of_Doubles.Split(xmin(i),x0,x1,x2,x3);
       x0 := x0 - m0;
       if x0 < 0.0
        then x0 := -x0;
       end if;
       x1 := x1 - m1;
       if x1 < 0.0
        then x1 := -x1;
       end if;
       x2 := x2 - m2;
       if x2 < 0.0
        then x2 := -x2;
       end if;
       x3 := x3 - m3;
       if x3 < 0.0
        then x3 := -x3;
       end if;
       xmbal := Balanced_Quarter_Doubles.Is_Balanced(0,x0,x1,x2,x3);
       put("x balanced : ");
       if xbal
        then put("true"); xcnt := xcnt + 1;
        else put("false");
       end if;
       put("  xplus balanced : ");
       if xpbal
        then put("true"); xpcnt := xpcnt + 1;
        else put("false");
       end if;
       put("  xmin balanced : ");
       if xmbal
        then put("true"); xmcnt := xpcnt + 1;
        else put("false");
       end if;
       new_line;
    end loop;
    put("    x balanced count : "); put(xcnt,1); new_line;
    put("xplus balanced count : "); put(xpcnt,1); new_line;
    put(" xmin balanced count : "); put(xmcnt,1); new_line;
  end Test_Balancing;

  procedure Main is

  -- DESCRIPTION :
  --   Launches the test on the inner product of doubles,
  --   balanced with a modifier constant.

    mct : constant double_float := Modifier_Constant(true);
    dim : integer32 := 0;

  begin
    put_line("Testing balanced representations of doubles ...");
    put("m : "); put(mct); new_line;
    put("Give the dimension : "); get(dim);
    Test_Inner_Product(dim);
    Test_Balancing(dim);
  end Main;

begin
  Main;
end ts_modbal;
