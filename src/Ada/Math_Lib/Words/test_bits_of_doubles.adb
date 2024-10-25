with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Natural_Vectors;
with Bits_of_Doubles;                    use Bits_of_Doubles;

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
    put("error :"); put(abs(y-x),2); new_line;
    z := x;
    chop_last_bits(z,26,thebits,lastbits);
    put("    z :"); put(z); new_line;
    put("last bits :"); write_52bits(lastbits); new_line;
    valbits := value_52bits(lastbits);
    put("h : "); put(valbits,1,b=>16); new_line;
  end Eat_Last_Bits_of_Pi;

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
    put("r :"); put(abs(x-z)); new_line;
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

  procedure Main is
  begin
    Eat_Last_Bits_of_Pi;
    Add_First_Bits_of_Pi;
  end Main;

end Test_Bits_of_Doubles;
