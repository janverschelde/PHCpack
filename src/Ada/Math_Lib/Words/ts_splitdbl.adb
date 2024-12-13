with text_io;                            use text_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Basics;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Bits_of_Doubles;                    use Bits_of_Doubles;

procedure ts_splitdbl is

-- DESCRIPTION :
--   Tests the splitting of a 64-bit floating-point number
--   in four equal parts, for use in making the product.
--   Two issues matter:
--   (1) A decimal representation may not have a finite binary expansion.
--   (2) A bit equal split does not suffice to give accurate results.

  procedure Test_Bit_Split ( x : in double_float ) is

  -- DESCRIPTION :
  --   Splits the fraction of x in four equal parts,
  --   using the bit representation of the fraction.

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
  end Test_Bit_Split;

  procedure to_Double_Double ( s : in string; x : out double_double;
                               verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Reads the string into a double double number
  --   to test for representation errors.

    fail : boolean;

  begin
    if verbose
     then put_line("parsing " & s & " ...");
    end if;
    Double_Double_Numbers_io.read(s,x,fail);
    if verbose then
      put("--> x : "); put(x);
      if fail
       then put_line(" FAILED!");
       else put_line(" okay");
      end if;
      put(" x hi :"); put(hi_part(x)); new_line;
      put(" x lo :"); put(lo_part(x)); new_line;
    end if;
  end to_Double_Double;

  procedure Test_Particular_Representations is

  -- DESCRIPTION :
  --   Tests the representations of two numbers,
  --   defined in decimal representation of double floats,
  --   which result in double doubles with a nonzero low double,
  --   because the decimal representation leads to a representation
  --   error when converted to binary.
  --   For the particular numbers in this test, multiplying only the
  --   high doubles leads to an error equal to machine precision,
  --   a representation error caused by omitting the low doubles.

    sx : constant string := "2.72666529413623016";
    sy : constant string := "2.09921004814059065000E-01";
    x,y,z : double_double;
    s32 : constant string := "5.7238431833669934787477014667440e-01";

  begin
    to_Double_Double(sx,x);
    to_Double_Double(sy,y);
    z := x*y;
    put("x*y : "); put(z); new_line;
    put("s32 : "); put(s32); new_line;
  end Test_Particular_Representations;

  procedure Test_Product ( x,y : in double_double ) is

  -- DESCRIPTION :
  --   Tests the product of two double doubles via splitting.

    z : constant double_double := x*y;
    p : double_double;
    x0,x1,x2,x3,x4,x5,x6,x7 : double_float;
    y0,y1,y2,y3,y4,y5,y6,y7 : double_float;
    z0,z1,z2,z3,z4,z5,z6,z7 : double_float;
    s0,s1,s2,s3,e1,e2,e3,e4 : double_float;

  begin
    Test_Bit_Split(hi_part(x));
    Test_Bit_Split(lo_part(x));
    Test_Bit_Split(hi_part(y));
    Test_Bit_Split(lo_part(y));
    Split(hi_part(x),x0,x1,x2,x3);
    Split(lo_part(x),x4,x5,x6,x7);
    Split(hi_part(y),y0,y1,y2,y3);
    Split(lo_part(y),y4,y5,y6,y7);
    put("x0 : "); write_52bits_expo(x0); new_line;
    put("x1 : "); write_52bits_expo(x1); new_line;
    put("x2 : "); write_52bits_expo(x2); new_line;
    put("x3 : "); write_52bits_expo(x3); new_line;
    put("x4 : "); write_52bits_expo(x4); new_line;
    put("x5 : "); write_52bits_expo(x5); new_line;
    put("x6 : "); write_52bits_expo(x6); new_line;
    put("x7 : "); write_52bits_expo(x7); new_line;
    put("y0 : "); write_52bits_expo(y0); new_line;
    put("y1 : "); write_52bits_expo(y1); new_line;
    put("y2 : "); write_52bits_expo(y2); new_line;
    put("y3 : "); write_52bits_expo(y3); new_line;
    put("y4 : "); write_52bits_expo(y4); new_line;
    put("y5 : "); write_52bits_expo(y5); new_line;
    put("y6 : "); write_52bits_expo(y6); new_line;
    put("y7 : "); write_52bits_expo(y7); new_line;
    z0 := x0*y0;
    z1 := x0*y1 + x1*y0;
    z2 := x0*y2 + x1*y1 + x2*y0;
    z3 := x0*y3 + x1*y2 + x2*y1 + x3*y0;
    z4 := x0*y4 + x1*y3 + x2*y2 + x3*y1 + x4*y0;
    z5 := x0*y5 + x1*y4 + x2*y3 + x3*y2 + x4*y1 + x5*y0;
    z6 := x0*y6 + x1*y5 + x2*y4 + x3*y3 + x4*y2 + x5*y1 + x6*y0;
    z7 := x0*y7 + x1*y6 + x2*y5 + x3*y4 + x4*y3 + x5*y2 + x6*y1 + x7*y0;
    put_line("The bits of the components of x*y :");
    put("z0 : "); write_52bits_expo(z0); new_line;
    put("z1 : "); write_52bits_expo(z1); new_line;
    put("z2 : "); write_52bits_expo(z2); new_line;
    put("z3 : "); write_52bits_expo(z3); new_line;
    put("z4 : "); write_52bits_expo(z4); new_line;
    put("z5 : "); write_52bits_expo(z5); new_line;
    put("z6 : "); write_52bits_expo(z6); new_line;
    put("z7 : "); write_52bits_expo(z7); new_line;
    Double_Double_Basics.two_sum(z0,z1,s0,e1);
    Double_Double_Basics.two_sum(z2,z3,s1,e2);
    Double_Double_Basics.two_sum(z4,z5,s2,e3);
    Double_Double_Basics.two_sum(z6,z7,s3,e4);
    put_line("after paired adding up :");
    put("s0 : "); write_52bits_expo(s0); new_line;
    put("e1 : "); write_52bits_expo(e1); new_line;
    put("s1 : "); write_52bits_expo(s1); new_line;
    put("e2 : "); write_52bits_expo(e2); new_line;
    put("s2 : "); write_52bits_expo(s2); new_line;
    put("e3 : "); write_52bits_expo(e3); new_line;
    put("s3 : "); write_52bits_expo(s3); new_line;
    put("e4 : "); write_52bits_expo(e4); new_line;
    p := Double_Double_Numbers.create(s3);
    p := p + s2;
    p := p + s1;
    p := p + s0;
    put("x*y : "); put(z); new_line;
    put("prd : "); put(p); new_line;
    put("phi : "); write_52bits_expo(hi_part(p)); new_line;
    put("zhi : "); write_52bits_expo(hi_part(z)); new_line;
    put("plo : "); write_52bits_expo(lo_part(p)); new_line;
    put("zlo : "); write_52bits_expo(lo_part(z)); new_line;
  end Test_Product;

  procedure Test_Particular_Product is

  -- DESCRIPTION :
  --   Tests the product of two particular doubles.

    sx : constant string := "2.72666529413623016";
    sy : constant string := "2.09921004814059065000E-01";
    x,y : double_double;
    s32 : constant string := "5.7238431833669934787477014667440e-01";

  begin
    to_Double_Double(sx,x);
    to_Double_Double(sy,y);
    Test_Product(x,y);
    put("s32 : "); put_line(s32);
  end Test_Particular_Product;

  procedure Main is

  -- DESCRIPTION :
  --   Runs tests on the split and the product of some particular numbers.

    x : constant double_float := 2.72666529413623016000E+00;

  begin
    Test_Bit_Split(x);
    Test_Particular_Representations;
    Test_Particular_Product;
  end Main;

begin
  Main;
end ts_splitdbl;
