with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Vectors;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Double_Double_Vectors;
with Double_Double_Vectors_io;           use Double_Double_Vectors_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Random_Vectors;
with Vectored_Double_Doubles;            use Vectored_Double_Doubles;

package body Test_Vectored_Double_Doubles is

  procedure Test_Real_Sum ( dim : in integer32 ) is

    ddz : constant DoblDobl_Complex_Vectors.Vector(1..dim)
        := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    ddv : Double_Double_Vectors.Vector(1..dim);
    v0,v1,v2,v3 : Standard_Floating_Vectors.Vector(1..dim);
    s0,s1,s2,s3 : double_float;
    ddsum0,ddsum1,err : double_double;

  begin
    for i in 1..dim loop
      ddv(i) := DoblDobl_Complex_Numbers.REAL_PART(ddz(i));
    end loop;
    put("Testing the sum of a random real vector of dimension ");
    put(dim,1);
    if dim > 20
     then put_line(" ...");
     else put_line(" :"); put_line(ddv);
    end if;
    Split(ddv,v0,v1,v2,v3);
    Sum(v0,v1,v2,v3,s0,s1,s2,s3);
    ddsum0 := Double_Double_Vectors.Sum(ddv);
    ddsum1 := Vectored_Double_Doubles.to_double_double(s0,s1,s2,s3);
    put("dd sum : "); put(ddsum0); new_line;
    put("dd vec : "); put(ddsum1); new_line;
    err := abs(ddsum0-ddsum1);
    put(" error : "); put(err,2); new_line;
  end Test_Real_Sum;

  procedure Test_Complex_Sum ( dim : in integer32 ) is

    ddv : constant DoblDobl_Complex_Vectors.Vector(1..dim)
        := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    v0re,v1re,v2re,v3re : Standard_Floating_Vectors.Vector(1..dim);
    v0im,v1im,v2im,v3im : Standard_Floating_Vectors.Vector(1..dim);
    s0re,s1re,s2re,s3re,s0im,s1im,s2im,s3im : double_float;
    ddsum0,ddsum1,err : Dobldobl_Complex_Numbers.Complex_Number;

    use DoblDobl_Complex_Numbers;

  begin
    put("Testing the sum of a random complex vector of dimension ");
    put(dim,1);
    if dim > 20
     then put_line(" ...");
     else put_line(" :"); put_line(ddv);
    end if;
    Split(ddv,v0re,v1re,v2re,v3re,v0im,v1im,v2im,v3im);
    Sum(v0re,v1re,v2re,v3re,v0im,v1im,v2im,v3im,
        s0re,s1re,s2re,s3re,s0im,s1im,s2im,s3im);
    ddsum0 := DoblDobl_Complex_Vectors.Sum(ddv);
    ddsum1 := to_Complex_Double_Double
                (s0re,s1re,s2re,s3re,s0im,s1im,s2im,s3im);
    put("dd sum : "); put(ddsum0); new_line;
    put("dd vec : "); put(ddsum1); new_line;
    err := ddsum0 - ddsum1;
    put(" error : "); put(err,2); new_line;
  end Test_Complex_Sum;

  procedure Test_Real_Product ( dim : in integer32 ) is

    a : constant DoblDobl_Complex_Vectors.Vector(1..dim)
      := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    b : constant DoblDobl_Complex_Vectors.Vector(1..dim)
      := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    x,y : Double_Double_Vectors.Vector(1..dim);
    x0,x1,x2,x3 : Standard_Floating_Vectors.Vector(1..dim);
    x4,x5,x6,x7 : Standard_Floating_Vectors.Vector(1..dim);
    y0,y1,y2,y3 : Standard_Floating_Vectors.Vector(1..dim);
    y4,y5,y6,y7 : Standard_Floating_Vectors.Vector(1..dim);
    s0,s1,s2,s3,s4,s5,s6,s7 : double_float;
    ddsum0,ddsum1,err : double_double;

  begin
    for i in 1..dim loop
      x(i) := DoblDobl_Complex_Numbers.REAL_PART(a(i));
      y(i) := DoblDobl_Complex_Numbers.REAL_PART(b(i));
    end loop;
    put("Testing product of random real vectors of dimension ");
    put(dim,1);
    if dim > 20 then
      put_line(" ...");
    else
      put_line(", x :"); put_line(x);
      put_line("y :"); put_line(y);
    end if;
    Quarter(x,x0,x1,x2,x3,x4,x5,x6,x7);
    Quarter(y,y0,y1,y2,y3,y4,y5,y6,y7);
    Product(x0,x1,x2,x3,x4,x5,x6,x7,y0,y1,y2,y3,y4,y5,y6,y7,
            s0,s1,s2,s3,s4,s5,s6,s7);
    ddsum0 := create(0.0);
    for i in x'range loop
      ddsum0 := ddsum0 + x(i)*y(i);
    end loop;
    ddsum1 := Vectored_Double_Doubles.to_double_double
                (s0,s1,s2,s3,s4,s5,s6,s7);
    put("dd prd : "); put(ddsum0); new_line;
    put("dd vec : "); put(ddsum1); new_line;
    err := abs(ddsum0-ddsum1);
    put(" error : "); put(err,2); new_line;
  end Test_Real_Product;

  procedure Test_Complex_Product ( dim : in integer32 ) is

    x : constant DoblDobl_Complex_Vectors.Vector(1..dim)
      := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    y : constant DoblDobl_Complex_Vectors.Vector(1..dim)
      := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    x0re,x1re,x2re,x3re : Standard_Floating_Vectors.Vector(1..dim);
    x4re,x5re,x6re,x7re : Standard_Floating_Vectors.Vector(1..dim);
    x0im,x1im,x2im,x3im : Standard_Floating_Vectors.Vector(1..dim);
    x4im,x5im,x6im,x7im : Standard_Floating_Vectors.Vector(1..dim);
    y0re,y1re,y2re,y3re : Standard_Floating_Vectors.Vector(1..dim);
    y4re,y5re,y6re,y7re : Standard_Floating_Vectors.Vector(1..dim);
    y0im,y1im,y2im,y3im : Standard_Floating_Vectors.Vector(1..dim);
    y4im,y5im,y6im,y7im : Standard_Floating_Vectors.Vector(1..dim);
    s0re,s1re,s2re,s3re,s4re,s5re,s6re,s7re : double_float;
    s0im,s1im,s2im,s3im,s4im,s5im,s6im,s7im : double_float;
    ddsum0,ddsum1,err : Dobldobl_Complex_Numbers.Complex_Number;

    use DoblDobl_Complex_Numbers;

  begin
    put("Testing on random complex vectors of dimension "); put(dim,1);
    if dim > 20 then
      put_line(" ...");
    else
      put_line(", x :"); put_line(x);
      put_line("y :"); put_line(y);
    end if;
    Quarter(x,x0re,x1re,x2re,x3re,x4re,x5re,x6re,x7re,
              x0im,x1im,x2im,x3im,x4im,x5im,x6im,x7im);
    Quarter(y,y0re,y1re,y2re,y3re,y4re,y5re,y6re,y7re,
              y0im,y1im,y2im,y3im,y4im,y5im,y6im,y7im);
    Product(x0re,x1re,x2re,x3re,x4re,x5re,x6re,x7re,
            x0im,x1im,x2im,x3im,x4im,x5im,x6im,x7im,
            y0re,y1re,y2re,y3re,y4re,y5re,y6re,y7re,
            y0im,y1im,y2im,y3im,y4im,y5im,y6im,y7im,
            s0re,s1re,s2re,s3re,s4re,s5re,s6re,s7re,
            s0im,s1im,s2im,s3im,s4im,s5im,s6im,s7im);
    ddsum0 := create(integer32(0));
    for i in x'range loop
      ddsum0 := ddsum0 + x(i)*y(i);
    end loop;
    ddsum1 := to_Complex_Double_Double
                (s0re,s1re,s2re,s3re,s4re,s5re,s6re,s7re,
                 s0im,s1im,s2im,s3im,s4im,s5im,s6im,s7im);
    put("dd prd : "); put(ddsum0); new_line;
    put("dd vec : "); put(ddsum1); new_line;
    err := ddsum0 - ddsum1;
    put(" error : "); put(err,2); new_line;
  end Test_Complex_Product;

  procedure Test_Complex_Norm ( dim : in integer32 ) is

    x : constant DoblDobl_Complex_Vectors.Vector(1..dim)
      := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    y : DoblDobl_Complex_Vectors.Vector(1..dim);
    x0re,x1re,x2re,x3re : Standard_Floating_Vectors.Vector(1..dim);
    x4re,x5re,x6re,x7re : Standard_Floating_Vectors.Vector(1..dim);
    x0im,x1im,x2im,x3im : Standard_Floating_Vectors.Vector(1..dim);
    x4im,x5im,x6im,x7im : Standard_Floating_Vectors.Vector(1..dim);
    y0re,y1re,y2re,y3re : Standard_Floating_Vectors.Vector(1..dim);
    y4re,y5re,y6re,y7re : Standard_Floating_Vectors.Vector(1..dim);
    y0im,y1im,y2im,y3im : Standard_Floating_Vectors.Vector(1..dim);
    y4im,y5im,y6im,y7im : Standard_Floating_Vectors.Vector(1..dim);
    s0re,s1re,s2re,s3re,s4re,s5re,s6re,s7re : double_float;
    s0im,s1im,s2im,s3im,s4im,s5im,s6im,s7im : double_float;
    ddsum0,ddsum1,err : Dobldobl_Complex_Numbers.Complex_Number;

    use DoblDobl_Complex_Numbers;

  begin
    for i in 1..dim loop
      y(i) := Conjugate(x(i));
    end loop;
    put("Testing on random complex vectors of dimension "); put(dim,1);
    if dim > 20 then
      put_line(" ...");
    else
      put_line(", x :"); put_line(x);
      put_line("y :"); put_line(y);
    end if;
    for i in 1..dim loop
      y(i) := Conjugate(x(i));
    end loop;
    Quarter(x,x0re,x1re,x2re,x3re,x4re,x5re,x6re,x7re,
              x0im,x1im,x2im,x3im,x4im,x5im,x6im,x7im);
    Quarter(y,y0re,y1re,y2re,y3re,y4re,y5re,y6re,y7re,
              y0im,y1im,y2im,y3im,y4im,y5im,y6im,y7im);
    Product(x0re,x1re,x2re,x3re,x4re,x5re,x6re,x7re,
            x0im,x1im,x2im,x3im,x4im,x5im,x6im,x7im,
            y0re,y1re,y2re,y3re,y4re,y5re,y6re,y7re,
            y0im,y1im,y2im,y3im,y4im,y5im,y6im,y7im,
            s0re,s1re,s2re,s3re,s4re,s5re,s6re,s7re,
            s0im,s1im,s2im,s3im,s4im,s5im,s6im,s7im);
    ddsum0 := create(integer32(0));
    for i in x'range loop
      ddsum0 := ddsum0 + x(i)*y(i);
    end loop;
    ddsum1 := to_Complex_Double_Double
                (s0re,s1re,s2re,s3re,s4re,s5re,s6re,s7re,
                 s0im,s1im,s2im,s3im,s4im,s5im,s6im,s7im);
    put("dd prd : "); put(ddsum0); new_line;
    put("dd vec : "); put(ddsum1); new_line;
    err := ddsum0 - ddsum1;
    put(" error : "); put(err,2); new_line;
  end Test_Complex_Norm;
  
  procedure Main is

    dim : integer32 := 0;

  begin
    put("Give the dimension : "); get(dim);
    Test_Real_Sum(dim);
    Test_Complex_Sum(dim);
    Test_Real_Product(dim);
    Test_Complex_Product(dim);
    Test_Complex_Norm(dim);
  end Main;

end Test_Vectored_Double_Doubles;