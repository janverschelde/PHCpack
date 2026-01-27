with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Random_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Standard_Floating_Vectors;
with Quad_Double_Vectors;
with Quad_Double_Vectors_io;             use Quad_Double_Vectors_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Random_Vectors;
with Vectored_Quad_Doubles;
with Random_Balanced_Quarters;

package body Test_Vectored_Quad_Doubles is

  procedure Test_Real_Product ( dim : in integer32 ) is

    a : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    b : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    xr,yr,x,y : Quad_Double_Vectors.Vector(1..dim);
    cntxpos,cntypos : natural32 := 0;
    qdprd0,qdprd0p,qdprd0m : quad_double;
    qdprd1,err : quad_double;

  begin
    for i in 1..dim loop
      xr(i) := QuadDobl_Complex_Numbers.REAL_PART(a(i));
      if xr(i) >= 0.0
       then cntxpos := cntxpos + 1;
      end if;
      yr(i) := QuadDobl_Complex_Numbers.REAL_PART(b(i));
      if yr(i) >= 0.0
       then cntypos := cntypos + 1;
      end if;
    end loop;
    x := Vectored_Quad_Doubles.Sign_Balance(xr,verbose=>false);
    y := Vectored_Quad_Doubles.Sign_Balance(yr,verbose=>false);
    put("Testing product of random real vectors of dimension ");
    put(dim,1);
    if dim > 20 then
      put_line(" ...");
    else
      put_line(", x :"); put_line(x);
      put_line("y :"); put_line(y);
    end if;
    put("# nonnegative x : "); put(cntxpos,1);
    put(", # nonnegative y : "); put(cntypos,1); new_line;
    qdprd0p := create(0.0);
    for i in x'range loop
      if (x(i) > 0.0 and y(i) > 0.0) or (x(i) < 0.0 and y(i) < 0.0)
       then qdprd0p := qdprd0p + x(i)*y(i);
      end if;
    end loop;
    qdprd0m := create(0.0);
    for i in x'range loop
      if (x(i) > 0.0 and y(i) < 0.0) or (x(i) < 0.0 and y(i) > 0.0)
       then qdprd0m := qdprd0m + x(i)*y(i);
      end if;
    end loop;
    qdprd0 := qdprd0p + qdprd0m;
    if dim > 20
     then qdprd1 := Vectored_Quad_Doubles.Product(x,y,false);
     else qdprd1 := Vectored_Quad_Doubles.Product(x,y);
    end if;
    put("qd prd : "); put(qdprd0); new_line;
    put("qd sgn : "); put(qdprd1); new_line;
    err := qdprd0 - qdprd1;
    put(" error : "); put(err,2); new_line;
  end Test_Real_Product;

  procedure Test_Complex_Product ( dim : in integer32 ) is

    x : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    y : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    qdprd0,qdprd1,err : Complex_Number;

  begin
    put("Testing on random complex vectors of dimension "); put(dim,1);
    if dim > 20 then
      put_line(" ...");
    else
      put_line(", x :"); put_line(x);
      put_line("y :"); put_line(y);
    end if;
    qdprd0 := create(integer32(0));
    for i in x'range loop
      qdprd0 := qdprd0 + x(i)*y(i);
    end loop;
    if dim > 20
     then qdprd1 := Vectored_Quad_Doubles.Product(x,y,false);
     else qdprd1 := Vectored_Quad_Doubles.Product(x,y);
    end if;
    put("qd prd : "); put(qdprd0); new_line;
    put("qd sgn : "); put(qdprd1); new_line;
    err := qdprd0 - qdprd1;
    put(" error : "); put(err,2); new_line;
  end Test_Complex_Product;

  procedure Test_Balanced_Product ( dim : in integer32 ) is

    x0,x1,x2,x3,x4,x5,x6,x7 : Standard_Floating_Vectors.Vector(1..dim);
    x8,x9,xA,xB,xC,xD,xE,xF : Standard_Floating_Vectors.Vector(1..dim);
    y0,y1,y2,y3,y4,y5,y6,y7 : Standard_Floating_Vectors.Vector(1..dim);
    y8,y9,yA,yB,yC,yD,yE,yF : Standard_Floating_Vectors.Vector(1..dim);
    x,y : Quad_Double_Vectors.Vector(1..dim);
    qdprd0,qdprd1,err : quad_double;
    s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,sA,sB,sC,sD,sE,sF : double_float;
    timer0,timer1 : Timing_Widget;
    freq : natural32 := 0;

  begin
    put_line("Testing the balanced inner product ...");
    put("Give the frequency : "); get(freq);
    Random_Balanced_Quarters.Random
      (dim,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,xA,xB,xC,xD,xE,xF);
    Random_Balanced_Quarters.Random
      (dim,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,yA,yB,yC,yD,yE,yF);
    x := Random_Balanced_Quarters.Make_Quad_Doubles
           (x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,xA,xB,xC,xD,xE,xF);
    y := Random_Balanced_Quarters.Make_Quad_Doubles
           (y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,yA,yB,yC,yD,yE,yF);
    tstart(timer0);
    for i in 1..freq loop
      qdprd0 := create(integer32(0));
      for i in x'range loop
        qdprd0 := qdprd0 + x(i)*y(i);
      end loop;
    end loop;
    tstop(timer0);
    tstart(timer1);
    for i in 1..freq loop
      Vectored_Quad_Doubles.Balanced_Quarter_Product
        (dim,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,xA,xB,xC,xD,xE,xF,
             y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,yA,yB,yC,yD,yE,yF,
         s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,sA,sB,sC,sD,sE,sF);
      if freq = 1 then
        qdprd1 := Vectored_Quad_Doubles.to_quad_double
          (s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,sA,sB,sC,sD,sE,sF);
      else
        qdprd1 := Vectored_Quad_Doubles.to_quad_double
          (s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,sA,sB,sC,sD,sE,sF,verbose=>false);
      end if;
    end loop;
    tstop(timer1);
    put("qd prd : "); put(qdprd0); new_line;
    put("qd sgn : "); put(qdprd1); new_line;
    err := qdprd0 - qdprd1;
    put(" error : "); put(err,2); new_line;
    new_line;
    print_times(standard_output,timer0,"quad double inner product");
    new_line;
    print_times(standard_output,timer1,"vectored quad double product");
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
   -- Test_Real_Product(dim);
   -- new_line;
   -- Test_Complex_Product(dim);
   -- new_line;
    Test_Balanced_Product(dim);
    put("Seed used : "); put(Standard_Random_Numbers.Get_Seed,1); new_line;
  end Main;

end Test_Vectored_Quad_Doubles;
