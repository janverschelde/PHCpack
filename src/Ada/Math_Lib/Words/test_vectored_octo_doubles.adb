with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Random_Numbers;
with Octo_Double_Numbers;                use Octo_Double_Numbers;
with Octo_Double_Numbers_io;             use Octo_Double_Numbers_io;
with OctoDobl_Complex_Numbers;           use OctoDobl_Complex_Numbers;
with OctoDobl_Complex_Numbers_io;        use OctoDobl_Complex_Numbers_io;
with Standard_Floating_Vectors;
with Octo_Double_Vectors;
with Octo_Double_Vectors_io;             use Octo_Double_Vectors_io;
with OctoDobl_Complex_Vectors;
with OctoDobl_Complex_Vectors_io;        use OctoDobl_Complex_Vectors_io;
with OctoDobl_Random_Vectors;
with Vectored_Octo_Doubles;
with Balanced_Quarter_Doubles;

package body Test_Vectored_Octo_Doubles is

  procedure Test_Real_Product ( dim : in integer32 ) is

    a : constant OctoDobl_Complex_Vectors.Vector(1..dim)
      := OctoDobl_Random_Vectors.Random_Vector(1,dim);
    b : constant OctoDobl_Complex_Vectors.Vector(1..dim)
      := OctoDobl_Random_Vectors.Random_Vector(1,dim);
    x,y : Octo_Double_Vectors.Vector(1..dim);
    cntxpos,cntypos : natural32 := 0;
    odprd0,odprd0p,odprd0m : Octo_double;
    odprd1,err : Octo_double;

  begin
    for i in 1..dim loop
      x(i) := OctoDobl_Complex_Numbers.REAL_PART(a(i));
      if x(i) >= 0.0
       then cntxpos := cntxpos + 1;
      end if;
      y(i) := OctoDobl_Complex_Numbers.REAL_PART(b(i));
      if y(i) >= 0.0
       then cntypos := cntypos + 1;
      end if;
    end loop;
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
    odprd0p := create(0.0);
    for i in x'range loop
      if (x(i) > 0.0 and y(i) > 0.0) or (x(i) < 0.0 and y(i) < 0.0)
       then odprd0p := odprd0p + x(i)*y(i);
      end if;
    end loop;
    odprd0m := create(0.0);
    for i in x'range loop
      if (x(i) > 0.0 and y(i) < 0.0) or (x(i) < 0.0 and y(i) > 0.0)
       then odprd0m := odprd0m + x(i)*y(i);
      end if;
    end loop;
    odprd0 := odprd0p + odprd0m;
    if dim > 20
     then odprd1 := Vectored_Octo_Doubles.Product(x,y,false);
     else odprd1 := Vectored_Octo_Doubles.Product(x,y);
    end if;
    put("od prd : "); put(odprd0); new_line;
    put("od sgn : "); put(odprd1); new_line;
    err := odprd0 - odprd1;
    put(" error : "); put(err,2); new_line;
  end Test_Real_Product;

  procedure Test_Complex_Product ( dim : in integer32 ) is

    x : constant OctoDobl_Complex_Vectors.Vector(1..dim)
      := OctoDobl_Random_Vectors.Random_Vector(1,dim);
    y : constant OctoDobl_Complex_Vectors.Vector(1..dim)
      := OctoDobl_Random_Vectors.Random_Vector(1,dim);
    odprd0 : Complex_Number;
    odprd1,err : Complex_Number;

  begin
    put("Testing on random complex vectors of dimension "); put(dim,1);
    if dim > 20 then
      put_line(" ...");
    else
      put_line(", x :"); put_line(x);
      put_line("y :"); put_line(y);
    end if;
    odprd0 := create(integer32(0));
    for i in x'range loop
      odprd0 := odprd0 + x(i)*y(i);
    end loop;
    if dim > 20
     then odprd1 := Vectored_Octo_Doubles.Product(x,y,false);
     else odprd1 := Vectored_Octo_Doubles.Product(x,y);
    end if;
    put("od prd : "); put(odprd0); new_line;
    put("od sgn : "); put(odprd1); new_line;
    err := odprd0 - odprd1;
    put(" error : "); put(err,2); new_line;
  end Test_Complex_Product;

  procedure Test_Balanced_Product ( dim : in integer32 ) is

    x00,x01,x02,x03 : Standard_Floating_Vectors.Vector(1..dim);
    x04,x05,x06,x07 : Standard_Floating_Vectors.Vector(1..dim);
    x08,x09,x10,x11 : Standard_Floating_Vectors.Vector(1..dim);
    x12,x13,x14,x15 : Standard_Floating_Vectors.Vector(1..dim);
    x16,x17,x18,x19 : Standard_Floating_Vectors.Vector(1..dim);
    x20,x21,x22,x23 : Standard_Floating_Vectors.Vector(1..dim);
    x24,x25,x26,x27 : Standard_Floating_Vectors.Vector(1..dim);
    x28,x29,x30,x31 : Standard_Floating_Vectors.Vector(1..dim);
    y00,y01,y02,y03 : Standard_Floating_Vectors.Vector(1..dim);
    y04,y05,y06,y07 : Standard_Floating_Vectors.Vector(1..dim);
    y08,y09,y10,y11 : Standard_Floating_Vectors.Vector(1..dim);
    y12,y13,y14,y15 : Standard_Floating_Vectors.Vector(1..dim);
    y16,y17,y18,y19 : Standard_Floating_Vectors.Vector(1..dim);
    y20,y21,y22,y23 : Standard_Floating_Vectors.Vector(1..dim);
    y24,y25,y26,y27 : Standard_Floating_Vectors.Vector(1..dim);
    y28,y29,y30,y31 : Standard_Floating_Vectors.Vector(1..dim);
    s00,s01,s02,s03,s04,s05,s06,s07 : double_float;
    s08,s09,s10,s11,s12,s13,s14,s15 : double_float;
    s16,s17,s18,s19,s20,s21,s22,s23 : double_float;
    s24,s25,s26,s27,s28,s29,s30,s31 : double_float;
    x,y : Octo_Double_Vectors.Vector(1..dim);
    odprd0,odprd1,err : octo_double;

  begin
    Balanced_Quarter_Doubles.Random
      (dim,x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x10,x11,x12,x13,x14,x15,
           x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31);
    Balanced_Quarter_Doubles.Random
      (dim,y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y10,y11,y12,y13,y14,y15,
           y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31);
    x := Balanced_Quarter_Doubles.Make_Octo_Doubles
           (x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x10,x11,x12,x13,x14,x15,
            x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31);
    y := Balanced_Quarter_Doubles.Make_Octo_Doubles
           (y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y10,y11,y12,y13,y14,y15,
            y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31);
    odprd0 := create(integer32(0));
    for i in x'range loop
      odprd0 := odprd0 + x(i)*y(i);
    end loop;
    Vectored_Octo_Doubles.Balanced_Quarter_Product
      (dim,x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x10,x11,x12,x13,x14,x15,
           x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,
           y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y10,y11,y12,y13,y14,y15,
           y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,
       s00,s01,s02,s03,s04,s05,s06,s07,s08,s09,s10,s11,s12,s13,s14,s15,
       s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,s29,s30,s31);
    odprd1 := Vectored_Octo_Doubles.to_octo_double
      (s00,s01,s02,s03,s04,s05,s06,s07,s08,s09,s10,s11,s12,s13,s14,s15,
       s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,s29,s30,s31);
    put("od prd : "); put(odprd0); new_line;
    put("od sgn : "); put(odprd1); new_line;
    err := odprd0 - odprd1;
    put(" error : "); put(err,2); new_line;
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
    Test_Real_Product(dim);
    new_line;
    Test_Complex_Product(dim);
    new_line;
    Test_Balanced_Product(dim);
    put("Seed used : "); put(Standard_Random_Numbers.Get_Seed,1); new_line;
  end Main;

end Test_Vectored_Octo_Doubles;
