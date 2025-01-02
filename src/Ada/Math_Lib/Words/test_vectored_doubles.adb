with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Random_Numbers;
with Standard_Floating_Vectors;
with Balanced_Quarter_Doubles;
with Vectored_Doubles;

package body Test_Vectored_Doubles is

  procedure Test_Balanced_Product ( dim : in integer32 ) is

    x0,x1,x2,x3 : Standard_Floating_Vectors.Vector(1..dim);
    y0,y1,y2,y3 : Standard_Floating_Vectors.Vector(1..dim);
    s0,s1,s2,s3 : double_float;
    x,y : Standard_Floating_Vectors.Vector(1..dim);
    dprd0,dprd1,err : double_float;
    timer0,timer1 : Timing_Widget;
    freq : natural32 := 0;
  
  begin
    put_line("Testing the balanced product ...");
    put("Give the frequency : "); get(freq);
    Balanced_Quarter_Doubles.Random(dim,x0,x1,x2,x3);
    Balanced_Quarter_Doubles.Random(dim,y0,y1,y2,y3);
    x := Balanced_Quarter_Doubles.Make_Doubles(x0,x1,x2,x3);
    y := Balanced_Quarter_Doubles.Make_Doubles(y0,y1,y2,y3);
    tstart(timer0);
    for i in 1..freq loop
      dprd0 := 0.0;
      for i in x'range loop
        dprd0 := dprd0 + x(i)*y(i);
      end loop;
    end loop;
    tstop(timer0);
    tstart(timer1);
    for i in 1..freq loop
      Vectored_Doubles.Balanced_Quarter_Product
        (dim,x0,x1,x2,x3,y0,y1,y2,y3,s0,s1,s2,s3);
      if freq = 1
       then dprd1 := Vectored_Doubles.to_double(s0,s1,s2,s3);
       else dprd1 := Vectored_Doubles.to_double(s0,s1,s2,s3,verbose=>false);
      end if;
    end loop;
    tstop(timer1);
    new_line;
    put("d prd : "); put(dprd0); new_line;
    put("d sgn : "); put(dprd1); new_line;
    err := dprd0 - dprd1;
    put("error : "); put(err,2); new_line;
    new_line;
    print_times(standard_output,timer0,"double inner product");
    new_line;
    print_times(standard_output,timer1,"vectored double product");
  end Test_Balanced_Product;

  procedure Main is

    seed : natural32 := 0;
    dim : integer32 := 0;

  begin
    put("Give the seed (0 for none) : "); get(seed);
    if seed /= 0
     then Standard_Random_Numbers.Set_Seed(seed);
    end if;
    put("Give the dimension : "); get(dim);
    Test_Balanced_Product(dim);
    put("Seed used : "); put(Standard_Random_Numbers.Get_Seed,1); new_line;
  end Main;

end Test_Vectored_Doubles;
