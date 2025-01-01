with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Floating_Vectors;
with Standard_Random_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Octo_Double_Numbers;                use Octo_Double_Numbers;
with Octo_Double_Numbers_io;             use Octo_Double_Numbers_io;
with Hexa_Double_Numbers;                use Hexa_Double_Numbers;
with Hexa_Double_Numbers_io;             use Hexa_Double_Numbers_io;
with Bits_of_Doubles;
with Balanced_Quarter_Doubles;

package body Test_Balanced_Quarter_Doubles is

  procedure Test_Thirteen_Bits is

    nbr : integer64;

  begin
    put_line("Making thirteen random bits ...");
    nbr := Balanced_Quarter_Doubles.Thirteen_Random_Bits;
    put("n : "); put(nbr,1); new_line;
    put("b : "); put(nbr,1,b=>2); new_line;
  end Test_Thirteen_Bits;

  procedure Write_Quarters ( x0,x1,x2,x3 : in double_float ) is

    x,s1,s2,s3,err : double_float;

  begin
    put("x0 :"); put(x0); new_line;
    put("x1 :"); put(x1); new_line;
    put("x2 :"); put(x2); new_line;
    put("x3 :"); put(x3); new_line;
    x := ((x3 + x2) + x1) + x0;
    put("          x : "); put(x); new_line;
    put("         x0 : "); put(x0); err := abs(x - x0);
    put("  error :"); put(err,2); new_line;
    s1 := x0 + x1;
    put("          x : "); put(x); new_line;
    put("      x0+x1 : "); put(s1); err := abs(x - s1);
    put("  error :"); put(err,2); new_line;
    s2 := x0 + x1 + x2;
    put("          x : "); put(x); new_line;
    put("   x0+x1+x2 : "); put(s2); err := abs(x - s2);
    put("  error :"); put(err,2); new_line;
    s3 := x0 + x1 + x2 + x3;
    put("          x : "); put(x); new_line;
    put("x0+x1+x2+x3 : "); put(s3); err := abs(x - s3);
    put("  error :"); put(err,2); new_line;
    put(" b0 : "); Bits_of_Doubles.write_52bits_expo(x0); new_line;
    put(" b1 : "); Bits_of_Doubles.write_52bits_expo(x1); new_line;
    put(" b2 : "); Bits_of_Doubles.write_52bits_expo(x2); new_line;
    put(" b3 : "); Bits_of_Doubles.write_52bits_expo(x3); new_line;
  end Write_Quarters;

  procedure Test_Random_Quarters is

    x0,x1,x2,x3 : double_float;

  begin
    put_line("Making a random balanced quarter double ...");
    Balanced_Quarter_Doubles.Random(x0,x1,x2,x3);
    Write_Quarters(x0,x1,x2,x3);
  end Test_Random_Quarters;

  procedure Test_Random_Vectors is

    dim : constant integer32 := 4;
    x0 : Standard_Floating_Vectors.Vector(1..dim);
    x1 : Standard_Floating_Vectors.Vector(1..dim);
    x2 : Standard_Floating_Vectors.Vector(1..dim);
    x3 : Standard_Floating_Vectors.Vector(1..dim);

  begin
    put_line("Making vectors of random balanced quarter doubles ...");
    Balanced_Quarter_Doubles.Random(dim,x0,x1,x2,x3);
    for i in 1..dim loop
      Write_Quarters(x0(i),x1(i),x2(i),x3(i));
    end loop;
  end Test_Random_Vectors;

  procedure Test_Double_Wrapper is

    x : constant double_float := Balanced_Quarter_Doubles.Random;
    x0,x1,x2,x3 : double_float;

  begin
    put("x : "); put(x); new_line;
    Bits_of_Doubles.Split(x,x0,x1,x2,x3);
    Write_Quarters(x0,x1,x2,x3);
  end Test_Double_Wrapper;

  procedure Test_Double_Double_Wrapper is

    x : constant double_double := Balanced_Quarter_Doubles.Random;
    x0,x1,x2,x3,x4,x5,x6,x7 : double_float;

  begin
    put("x : "); put(x); new_line;
    Bits_of_Doubles.Split(hi_part(x),x0,x1,x2,x3);
    Bits_of_Doubles.Split(lo_part(x),x4,x5,x6,x7);
    Write_Quarters(x0,x1,x2,x3);
    Write_Quarters(x4,x5,x6,x7);
  end Test_Double_Double_Wrapper;

  procedure Test_Quad_Double_Wrapper is

    x : constant quad_double := Balanced_Quarter_Doubles.Random;
    x0,x1,x2,x3,x4,x5,x6,x7 : double_float;
    x8,x9,xA,xB,xC,xD,xE,xF : double_float;

  begin
    put("x : "); put(x); new_line;
    Bits_of_Doubles.Split(hihi_part(x),x0,x1,x2,x3);
    Bits_of_Doubles.Split(lohi_part(x),x4,x5,x6,x7);
    Bits_of_Doubles.Split(hilo_part(x),x8,x9,xA,xB);
    Bits_of_Doubles.Split(lolo_part(x),xC,xD,xE,xF);
    Write_Quarters(x0,x1,x2,x3);
    Write_Quarters(x4,x5,x6,x7);
    Write_Quarters(x8,x9,xA,xB);
    Write_Quarters(xC,xD,xE,xF);
  end Test_Quad_Double_Wrapper;

  procedure Test_Octo_Double_Wrapper is

    x : constant octo_double := Balanced_Quarter_Doubles.Random;
    x00,x01,x02,x03,x04,x05,x06,x07 : double_float;
    x08,x09,x10,x11,x12,x13,x14,x15 : double_float;
    x16,x17,x18,x19,x20,x21,x22,x23 : double_float;
    x24,x25,x26,x27,x28,x29,x30,x31 : double_float;

  begin
    put("x : "); put(x); new_line;
    Bits_of_Doubles.Split(hihihi_part(x),x00,x01,x02,x03);
    Bits_of_Doubles.Split(lohihi_part(x),x04,x05,x06,x07);
    Bits_of_Doubles.Split(hilohi_part(x),x08,x09,x10,x11);
    Bits_of_Doubles.Split(lolohi_part(x),x12,x13,x14,x15);
    Bits_of_Doubles.Split(hihilo_part(x),x16,x17,x18,x19);
    Bits_of_Doubles.Split(lohilo_part(x),x20,x21,x22,x23);
    Bits_of_Doubles.Split(hilolo_part(x),x24,x25,x26,x27);
    Bits_of_Doubles.Split(lololo_part(x),x28,x29,x30,x31);
    Write_Quarters(x00,x01,x02,x03);
    Write_Quarters(x04,x05,x06,x07);
    Write_Quarters(x08,x09,x10,x11);
    Write_Quarters(x12,x13,x14,x15);
    Write_Quarters(x16,x17,x18,x19);
    Write_Quarters(x20,x21,x22,x23);
    Write_Quarters(x24,x25,x26,x27);
    Write_Quarters(x28,x29,x30,x31);
  end Test_Octo_Double_Wrapper;

  procedure Test_Hexa_Double_Wrapper is

    x : constant hexa_double := Balanced_Quarter_Doubles.Random;
    x00,x01,x02,x03,x04,x05,x06,x07 : double_float;
    x08,x09,x10,x11,x12,x13,x14,x15 : double_float;
    x16,x17,x18,x19,x20,x21,x22,x23 : double_float;
    x24,x25,x26,x27,x28,x29,x30,x31 : double_float;
    x32,x33,x34,x35,x36,x37,x38,x39 : double_float;
    x40,x41,x42,x43,x44,x45,x46,x47 : double_float;
    x48,x49,x50,x51,x52,x53,x54,x55 : double_float;
    x56,x57,x58,x59,x60,x61,x62,x63 : double_float;

  begin
    put("x : "); put(x); new_line;
    Bits_of_Doubles.Split(hihihihi_part(x),x00,x01,x02,x03);
    Bits_of_Doubles.Split(lohihihi_part(x),x04,x05,x06,x07);
    Bits_of_Doubles.Split(hilohihi_part(x),x08,x09,x10,x11);
    Bits_of_Doubles.Split(lolohihi_part(x),x12,x13,x14,x15);
    Bits_of_Doubles.Split(hihilohi_part(x),x16,x17,x18,x19);
    Bits_of_Doubles.Split(lohilohi_part(x),x20,x21,x22,x23);
    Bits_of_Doubles.Split(hilolohi_part(x),x24,x25,x26,x27);
    Bits_of_Doubles.Split(lololohi_part(x),x28,x29,x30,x31);
    Bits_of_Doubles.Split(hihihilo_part(x),x32,x33,x34,x35);
    Bits_of_Doubles.Split(lohihilo_part(x),x36,x37,x38,x39);
    Bits_of_Doubles.Split(hilohilo_part(x),x40,x41,x42,x43);
    Bits_of_Doubles.Split(lolohilo_part(x),x44,x45,x46,x47);
    Bits_of_Doubles.Split(hihilolo_part(x),x48,x49,x50,x51);
    Bits_of_Doubles.Split(lohilolo_part(x),x52,x53,x54,x55);
    Bits_of_Doubles.Split(hilololo_part(x),x56,x57,x58,x59);
    Bits_of_Doubles.Split(lolololo_part(x),x60,x61,x62,x63);
    Write_Quarters(x00,x01,x02,x03);
    Write_Quarters(x04,x05,x06,x07);
    Write_Quarters(x08,x09,x10,x11);
    Write_Quarters(x12,x13,x14,x15);
    Write_Quarters(x16,x17,x18,x19);
    Write_Quarters(x20,x21,x22,x23);
    Write_Quarters(x24,x25,x26,x27);
    Write_Quarters(x28,x29,x30,x31);
    Write_Quarters(x32,x33,x34,x35);
    Write_Quarters(x36,x37,x38,x39);
    Write_Quarters(x40,x41,x42,x43);
    Write_Quarters(x44,x45,x46,x47);
    Write_Quarters(x48,x49,x50,x51);
    Write_Quarters(x52,x53,x54,x55);
    Write_Quarters(x56,x57,x58,x59);
    Write_Quarters(x60,x61,x62,x63);
  end Test_Hexa_Double_Wrapper;

  procedure Test_Balanced_Split is

    x : constant double_float := Standard_Random_Numbers.Random;
    x0,x1,x2,x3 : double_float;
    isbal : boolean;

  begin
    put("x : "); put(x); new_line;
    put("b : "); Bits_of_Doubles.write_52bits_expo(x); new_line;
    Bits_of_Doubles.Split(x,x0,x1,x2,x3);
    Write_Quarters(x0,x1,x2,x3);
    isbal := Balanced_Quarter_Doubles.Is_Balanced(0,x0,x1,x2,x3);
    if not isbal
     then put_line("unbalanced");
     else put_line("balanced");
    end if;
  end Test_Balanced_Split;

  procedure Main is
  begin
    Test_Thirteen_Bits;
    Test_Random_Quarters;
    Test_Random_Vectors;
    new_line;
    put_line("Testing double wrapper ...");
    Test_Double_Wrapper;
    new_line;
    put_line("Testing double double wrapper ...");
    Test_Double_Double_Wrapper;
    new_line;
    put_line("Testing quad double wrapper ...");
    Test_Quad_Double_Wrapper;
    new_line;
    put_line("Testing octo double wrapper ...");
    Test_Octo_Double_Wrapper;
    new_line;
    put_line("Testing hexa double wrapper ...");
    Test_Hexa_Double_Wrapper;
    new_line;
    put_line("Testing the balanced split ...");
    Test_Balanced_Split;
  end Main;

end Test_Balanced_Quarter_Doubles;
