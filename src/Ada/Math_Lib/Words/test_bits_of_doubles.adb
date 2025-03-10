with text_io;                            use text_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Floating_Vectors;
with Standard_Random_Vectors;
with Double_Double_Basics;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Bits_of_Doubles;                    use Bits_of_Doubles;
with Mask_Bits_of_Doubles;

package body Test_Bits_of_Doubles is

  procedure Eat_Last_Bits_of_Pi is

    x : constant double_float := 3.141592653589793;
    f : constant double_float := double_float'fraction(x);
    e : constant integer32 := integer32(double_float'exponent(x));
    c : constant double_float := double_float'compose(f, e);
    s : constant double_float := double_float'compose(f, 52);
    m : constant integer64 := integer64(double_float'truncation(s));
    thebits,lastbits : Standard_Natural_Vectors.Vector(0..51);
    valbits : integer64;
    y,z : double_float;

  begin
    new_line;
    put_line("*** removing the last 26 bits of pi ***");
    new_line;
    put("  The 64-bit float : ");
    put("x : "); put(x); new_line;
    put("->    the fraction : ");
    put("f : "); put(f); new_line;
    put("->    the exponent :");
    put(" e : "); put(e); new_line;
    put("-> composed number : ");
    put("c : "); put(c); new_line;
    put("-> 52-bit fraction : ");
    put("m : "); put(m); new_line;
    expand_52bits(thebits,m);
    put(" the bits :"); write_52bits(thebits); new_line;
    valbits := value_52bits(thebits);
    put("->    52-bit value : ");
    put("v : "); put(valbits); new_line;
    put("h : "); put(valbits,1,b=>16); new_line;
    put("b : "); put(valbits,1,b=>2); new_line;
    for i in 0..25 loop
      thebits(integer32(51-i)) := 0;
    end loop;
    put(" the bits :"); write_52bits(thebits); new_line;
    valbits := value_52bits(thebits);
    put("->    52-bit value : ");
    put("v : "); put(valbits); new_line;
    put("h : "); put(valbits,1,b=>16); new_line;
    put("b : "); put(valbits,1,b=>2); new_line;
    y := double_float'compose(double_float(valbits), e);
    z := chop_last_bits(x,26);
    put("    x :"); put(x); new_line;
    put("    y :"); put(y); new_line;
    put("    z :"); put(z); new_line;
    put("error :"); Standard_Floating_Numbers_io.put(abs(y-x),2); new_line;
    z := x;
    chop_last_bits(z,26,thebits,lastbits);
    put("    z :"); put(z); new_line;
    put("last bits :"); write_52bits(lastbits); new_line;
    valbits := value_52bits(lastbits);
    put("h : "); put(valbits,1,b=>16); new_line;
  end Eat_Last_Bits_of_Pi;

  procedure Mod_Last_Bits_of_Pi is

    x : constant double_float := 3.141592653589793;
    f : constant double_float := double_float'fraction(x);
    e : constant integer32 := integer32(double_float'exponent(x));
    c : constant double_float := double_float'compose(f, e);
    s : constant double_float := double_float'compose(f, 52);
    m : constant integer64 := integer64(double_float'truncation(s));
    mlast,mchop : integer64;
    y : double_float;

  begin
    new_line;
    put_line("*** modding out the last 26 bits of pi ***");
    new_line;
    put("  The 64-bit float : ");
    put("x : "); put(x); new_line;
    put("->    the fraction : ");
    put("f : "); put(f); new_line;
    put("->    the exponent :");
    put(" e : "); put(e); new_line;
    put("-> composed number : ");
    put("c : "); put(c); new_line;
    put("-> 52-bit fraction : ");
    put("m : "); put(m); new_line;
    put_line(" the bits :"); 
    put("m : "); put(m,1,b=>2); new_line;
    put("m : "); put(m,1,b=>16); new_line;
    mlast := Mask_Bits_of_Doubles.last_bits(m,26);
    put("last bits : "); put(mlast,1,b=>2); new_line;
    put("last hex  : "); put(mlast,1,b=>16); new_line;
    mchop := m - mlast;
    put("chop bits : "); put(mchop,1,b=>2); new_line;
    put("chop hex  : "); put(mchop,1,b=>16); new_line;
    y := double_float'compose(double_float(mchop),e);
    put("    x :"); put(x); new_line;
    put("    y :"); put(y); new_line;
    put("error :"); Standard_Floating_Numbers_io.put(abs(y-x),2); new_line;
  end Mod_Last_Bits_of_Pi;

  procedure Add_First_Bits_of_Pi is

    x : constant double_float := 3.141592653589793;
    e : constant integer32 := integer32(double_float'exponent(x));
    firstbits,lastbits : Standard_Natural_Vectors.Vector(0..51);
    thebits : Standard_Natural_Vectors.Vector(0..51);
    valbits : integer64;
    z : double_float := x;
    y,fy,sy,w : double_float;
    ey : integer32;
    my : integer64;

  begin
    new_line;
    put_line("*** adding the first 26 bits to the last 26 bits of pi ***");
    new_line;
    chop_last_bits(z,26,firstbits,lastbits);
    valbits := value_52bits(lastbits);
    y := double_float'compose(double_float(valbits), e-26);
    put("  The 64-bit float : ");
    put("x :"); put(x); new_line;
    put(" The first 26 bits : ");
    put("z :"); put(z); new_line;
    put("  The last 26 bits : ");
    put("y :"); put(y); new_line;
    put("         The error : ");
    put("r :"); Standard_Floating_Numbers_io.put(abs(x-z)); new_line;
    fy := double_float'fraction(y);
    ey := integer32(double_float'exponent(y));
    put("e : "); put(ey); new_line;
    sy := double_float'compose(fy, 52);
    my := integer64(double_float'truncation(sy));
    put("m : "); put(my,1,b=>2); new_line;
    put(" last bits :"); write_52bits(lastbits); new_line;
    put("first bits :"); write_52bits(firstbits); new_line;
    expand_52bits(thebits,my);
    put("  the bits :"); write_52bits(thebits); new_line;
    for i in 0..25 loop
      thebits(integer32(51-i)) := thebits(integer32(51-i-26));
      thebits(integer32(51-i-26)) := 0;
    end loop;
    put("   shifted :"); write_52bits(thebits); new_line;
    for i in 0..25 loop
      thebits(integer32(i)) := firstbits(integer32(i));
    end loop;
    put(" add first :"); write_52bits(thebits); new_line;
    valbits := value_52bits(thebits);
    z := double_float'compose(double_float(valbits),ey + 26);
    put("z :"); put(z); new_line;
    put("x :"); put(x); new_line;
    put("y :"); put(y); new_line;
    w := Insert_First_Bits(y,26,firstbits);
    put("w :"); put(w); new_line;
    put("x :"); put(x); new_line;
  end Add_First_Bits_of_Pi;

  procedure Mod_First_Bits_of_Pi is

    x : constant double_float := 3.141592653589793;
    e : constant integer32 := integer32(double_float'exponent(x));
    f : constant double_float := double_float'fraction(x);
    s : constant double_float := double_float'compose(f, 52);
    m : constant integer64 := integer64(double_float'truncation(s));
    mlast,mchop,mfirst : integer64;
    z : double_float := x;
    y,w,a : double_float;

  begin
    new_line;
    put_line("*** modding the first 26 bits to the last 26 bits of pi ***");
    new_line;
    mlast := Mask_Bits_of_Doubles.last_bits(m,26);
    put("last bits : "); put(mlast,1,b=>2); new_line;
    put("last hex  : "); put(mlast,1,b=>16); new_line;
    mchop := m - mlast;
    put("chop bits : "); put(mchop,1,b=>2); new_line;
    put("chop hex  : "); put(mchop,1,b=>16); new_line;
    y := double_float'compose(double_float(mlast), e-26);
    z := double_float'compose(double_float(mchop), e);
    put("  The 64-bit float : ");
    put("x :"); put(x); new_line;
    put(" The first 26 bits : ");
    put("z :"); put(z); new_line;
    put("  The last 26 bits : ");
    put("y :"); put(y); new_line;
    put("         The error : ");
    put("r :"); Standard_Floating_Numbers_io.put(abs(x-z)); new_line;
    mfirst := Mask_Bits_of_Doubles.first_bits(m,26);
    put("first bits : "); put(mfirst,1,b=>2); new_line;
    put("first bits : "); put(mfirst,1,b=>16); new_line;
    w := double_float'compose(double_float(mfirst), e);
    put(" The first 26 bits : ");
    put("w :"); put(w); new_line;
    a := w + y;
    put("  Added to last 26 : ");
    put("a :"); put(a); new_line;
    put("  The 64-bit float : ");
    put("x :"); put(x); new_line;
  end Mod_First_Bits_of_Pi;

  procedure Test_Mod_Mask_Bits is
  begin
    Eat_Last_Bits_of_Pi;
    Mod_Last_Bits_of_Pi;
    Add_First_Bits_of_Pi;
    Mod_First_Bits_of_Pi;
  end Test_Mod_Mask_Bits;

  procedure Test_Last_Bits ( lastbits : in natural32 ) is

    use Mask_Bits_of_Doubles;

    r : constant unsigned_integer64 := 2**natural(lastbits);
    x : unsigned_integer64 := r;
    m : constant unsigned_integer64 := last_mask(lastbits);
    y : unsigned_integer64;

  begin
    put("   x : "); put(integer64(x),1,b=>2); new_line;
    put("mask :  "); put(integer64(m),1,b=>2); new_line;
    for i in 1..r loop
      y := x and m;
      put("last bits of "); put(integer64(x)); put(" : ");
      put(integer64(y),1,b=>2); put(" = ");
      put(last_bits(integer64(x),lastbits),1,b=>2);
      if y = last_bits(x,lastbits)
       then put_line(" true");
       else put_line(" false");
      end if;
      x := x + 1;
    end loop;
  end Test_Last_Bits;

  procedure Test_First_Bits ( firstbits : in natural32 ) is

    use Mask_Bits_of_Doubles;

    r : constant unsigned_integer64 := 2**natural(firstbits);
    xexp : constant integer32 := 51 - integer32(firstbits);
    xfloor : constant unsigned_integer64 := 2**natural(xexp);
    modulus : constant unsigned_integer64 := 2*xfloor;
    x : unsigned_integer64 := xfloor;
    y,z : unsigned_integer64;

  begin
    put("2^"); put(xexp); put(" : ");
    put(integer64(xfloor),1,b=>2); new_line;
    put("   x : "); put(integer64(x),1,b=>2); new_line;
    for i in 0..(r-1) loop
     -- y := (x - last_bits(x,52 - firstbits));
      y := first_bits(x,firstbits);
      z := y/modulus;
      put(integer64(x),59,b=>2); new_line;
      put(integer64(y),59,b=>2); new_line;
      put(integer64(z),10,b=>2);
      put(" = "); put(integer32(i));
      if z = i
       then put_line(" true");
       else put_line(" false");
      end if;
      x := x + modulus;
    end loop;
  end Test_First_bits;

  procedure Test_Mask_Bits is
  begin
    for i in 1..4 loop
      put("testing the extraction of the last ");
      put(integer32(i),1); put_line(" bits ...");
      test_last_bits(natural32(i));
    end loop;
    for i in 1..4 loop
      put("testing the extraction of the first ");
      put(integer32(i),1); put_line(" bits ...");
      test_first_bits(natural32(i));
    end loop;
  end Test_Mask_Bits;

  procedure Test_Bit_Split ( x : in double_float ) is

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

  procedure Test_Split_Product is

    x : constant double_float := 2.72666529413623016000E+00;

  begin
    Test_Bit_Split(x);
    Test_Particular_Representations;
    Test_Particular_Product;
  end Test_Split_Product;

  procedure Test_Last_Zero_Counts is

    x : constant double_float := abs(Standard_Random_Numbers.Random);
    nbr : natural32;
    x0,x1,x2,x3 : double_float;

  begin
    put("x :"); put(x); new_line;
    put("b : "); write_fraction_bits(x);
    nbr := Last_Zero_Count(x);
    put("n : "); put(nbr,1); new_line;
    put_line("after splitting in half : ");
    Mod_Split(x,x0,x1,verbose=>false);
    put("x0 : "); write_fraction_bits(x0);
    put("n0 : "); put(Last_Zero_Count(x0),1); new_line;
    put("x1 : "); write_fraction_bits(x1);
    put("n1 : "); put(Last_Zero_Count(x1),1); new_line;
    put_line("after quartering : ");
    Split(x,x0,x1,x2,x3);
    put("x0 : "); write_fraction_bits(x0);
    put("n0 : "); put(Last_Zero_Count(x0),1); new_line;
    put("x1 : "); write_fraction_bits(x1);
    put("n1 : "); put(Last_Zero_Count(x1),1); new_line;
    put("x2 : "); write_fraction_bits(x2);
    put("n2 : "); put(Last_Zero_Count(x2),1); new_line;
    put("x3 : "); write_fraction_bits(x3);
    put("n3 : "); put(Last_Zero_Count(x3),1); new_line;
  end Test_Last_Zero_Counts;

  procedure Test_Free_Convolution_Bits ( dim : in integer32 ) is

    x : Standard_Floating_Vectors.Vector(1..dim)
      := Standard_Random_Vectors.Random_Vector(1,dim);
    y : Standard_Floating_Vectors.Vector(1..dim)
      := Standard_Random_Vectors.Random_Vector(1,dim);
    x0,x1,x2,x3 : Standard_Floating_Vectors.Vector(1..dim);
    y0,y1,y2,y3 : Standard_Floating_Vectors.Vector(1..dim);
    s0,s1,s2,s3 : double_float;
    f : integer64;
    e : integer32;
    xemin0,xemin1,xemin2,xemin3 : integer32;
    xemax0,xemax1,xemax2,xemax3 : integer32;
    yemin0,yemin1,yemin2,yemin3 : integer32;
    yemax0,yemax1,yemax2,yemax3 : integer32;

    procedure Update ( f : in integer64;
                       emin,emax : in out integer32; e : in integer32 ) is

    -- DESCRIPTION :
    --   Updates the minimal and maximal exponent with e
    --   if the fraction is nonzero.

    begin
      if f /= 0 then
        if e < emin
         then emin := e;
        end if;
        if e > emax
         then emax := e;
        end if;
      end if;
    end Update;

  begin
    x(1) := abs(x(1)); Split(x(1),x0(1),x1(1),x2(1),x3(1));
    Fraction_Exponent(x0(1),f,xemin0); xemax0 := xemin0;
    Fraction_Exponent(x1(1),f,xemin1); xemax1 := xemin1;
    Fraction_Exponent(x2(1),f,xemin2); xemax2 := xemin2;
    Fraction_Exponent(x3(1),f,xemin3); xemax3 := xemin3;
    y(1) := abs(y(1)); Split(y(1),y0(1),y1(1),y2(1),y3(1));
    Fraction_Exponent(y0(1),f,yemin0); yemax0 := yemin0;
    Fraction_Exponent(y1(1),f,yemin1); yemax1 := yemin1;
    Fraction_Exponent(y2(1),f,yemin2); yemax2 := yemin2;
    Fraction_Exponent(y3(1),f,yemin3); yemax3 := yemin3;
    for i in 2..dim loop
      x(i) := abs(x(i)); Split(x(i),x0(i),x1(i),x2(i),x3(i));
      Fraction_Exponent(x0(i),f,e); Update(f,xemin0,xemax0,e);
      Fraction_Exponent(x1(i),f,e); Update(f,xemin1,xemax1,e);
      Fraction_Exponent(x2(i),f,e); Update(f,xemin2,xemax2,e);
      Fraction_Exponent(x3(i),f,e); Update(f,xemin3,xemax3,e);
      y(i) := abs(y(i)); Split(y(i),y0(i),y1(i),y2(i),y3(i));
      Fraction_Exponent(y0(i),f,e); Update(f,yemin0,yemax0,e);
      Fraction_Exponent(y1(i),f,e); Update(f,yemin1,yemax1,e);
      Fraction_Exponent(y2(i),f,e); Update(f,yemin2,yemax2,e);
      Fraction_Exponent(y3(i),f,e); Update(f,yemin3,yemax3,e);
    end loop;
    s0 := 0.0; s1 := 0.0; s2 := 0.0; s3 := 0.0;
    for i in 1..dim loop
      s0 := s0 + x0(i)*y0(i);
      s1 := s1 + x0(i)*y1(i) + x1(i)*y0(i);
      s2 := s2 + x0(i)*y2(i) + x1(i)*y1(i) + x2(i)*y0(i);
      s3 := s3 + x0(i)*y3(i) + x1(i)*y2(i) + x2(i)*y1(i) + x3(i)*y0(i);
    end loop;
    put("n0 : "); put(Last_Zero_Count(s0),1);
    put(", n1 : "); put(Last_Zero_Count(s1),1); 
    put(", n2 : "); put(Last_Zero_Count(s2),1);
    put(", n3 : "); put(Last_Zero_Count(s3),1); new_line;
    put("  xemin : ");
    put(" "); put(xemin0,1); 
    put(" "); put(xemin1,1); 
    put(" "); put(xemin2,1); 
    put(" "); put(xemin3,1);
    put(", xemax : ");
    put(" "); put(xemax0,1); 
    put(" "); put(xemax1,1); 
    put(" "); put(xemax2,1); 
    put(" "); put(xemax3,1);
    put(", d : "); put(xemax0-xemin0,1);
    put(" "); put(xemax1-xemin1,1);
    put(" "); put(xemax2-xemin2,1);
    put(" "); put(xemax3-xemin3,1); new_line;
    put("  yemin : ");
    put(" "); put(yemin0,1); 
    put(" "); put(yemin1,1); 
    put(" "); put(yemin2,1); 
    put(" "); put(yemin3,1);
    put(", yemax : ");
    put(" "); put(yemax0,1); 
    put(" "); put(yemax1,1); 
    put(" "); put(yemax2,1); 
    put(" "); put(yemax3,1);
    put(", d : "); put(yemax0-yemin0,1);
    put(" "); put(yemax1-yemin1,1);
    put(" "); put(yemax2-yemin2,1);
    put(" "); put(yemax3-yemin3,1); new_line;
  end Test_Free_Convolution_Bits;

  procedure Main is

    dim : integer32 := 0;

  begin
    Test_Mod_Mask_Bits;
    new_line;
    put_line("*** testing the masking of bits of doubles ***");
    new_line;
    Test_Mask_Bits;
    new_line;
    put_line("*** testing the product via splits ***");
    new_line;
    Test_Split_Product;
    new_line;
    put_line("*** testing the last zero counts ***");
    Test_Last_Zero_Counts;
    new_line;
    put("Give the dimension : "); get(dim);
    for i in 1..10 loop
      Test_Free_Convolution_Bits(dim);
    end loop;
  end Main;

end Test_Bits_of_Doubles;
