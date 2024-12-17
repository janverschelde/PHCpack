with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
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

procedure ts_vdda is

-- DESCRIPTION :
--   Test on the vectorized double double arithmetic.

  procedure Test_Real_Sum ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random vector of double double numbers
  --   of dimension dim.

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
    put("Testing on a random vector of dimension "); put(dim,1);
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

  -- DESCRIPTION :
  --   Generates a random vector of double double numbers
  --   of dimension dim.

    ddv : constant DoblDobl_Complex_Vectors.Vector(1..dim)
        := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    v0re,v1re,v2re,v3re : Standard_Floating_Vectors.Vector(1..dim);
    v0im,v1im,v2im,v3im : Standard_Floating_Vectors.Vector(1..dim);
    s0re,s1re,s2re,s3re,s0im,s1im,s2im,s3im : double_float;
    ddsum0,ddsum1,err : Dobldobl_Complex_Numbers.Complex_Number;

    use DoblDobl_Complex_Numbers;

  begin
    put("Testing on a random complex vector of dimension ");
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
  
  procedure Main is

  -- DESCRIPTION :
  --   Prompts the dimension and runs the tests.

    dim : integer32 := 0;

  begin
    put("Give the dimension : "); get(dim);
    Test_Real_Sum(dim);
    Test_Complex_Sum(dim);
  end Main;

begin
  Main;
end ts_vdda;
