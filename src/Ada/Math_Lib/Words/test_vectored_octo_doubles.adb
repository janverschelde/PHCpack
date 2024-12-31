with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Random_Numbers;
with Octo_Double_Numbers;                use Octo_Double_Numbers;
with Octo_Double_Numbers_io;             use Octo_Double_Numbers_io;
with OctoDobl_Complex_Numbers;           use OctoDobl_Complex_Numbers;
with OctoDobl_Complex_Numbers_io;        use OctoDobl_Complex_Numbers_io;
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

    x : constant Octo_Double_Vectors.Vector(1..dim)
      := Balanced_Quarter_Doubles.Random(dim);
    y : constant Octo_Double_Vectors.Vector(1..dim)
      := Balanced_Quarter_Doubles.Random(dim);
    odprd0,odprd1,err : octo_double;

  begin
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
