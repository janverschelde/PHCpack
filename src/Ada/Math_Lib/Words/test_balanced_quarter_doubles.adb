with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Floating_Vectors;
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

  procedure Main is
  begin
    Test_Thirteen_Bits;
    Test_Random_Quarters;
    Test_Random_Vectors;
  end Main;

end Test_Balanced_Quarter_Doubles;
