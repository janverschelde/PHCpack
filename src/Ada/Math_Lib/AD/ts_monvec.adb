with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Vector_Norms;
with Standard_Random_Vectors;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;       use DoblDobl_Complex_Numbers_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;       use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_Vector_Norms;
with DoblDobl_Random_Vectors;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;       use QuadDobl_Complex_Numbers_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;       use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vector_Norms;
with QuadDobl_Random_Vectors;
with Standard_Monomial_Vectors;
with Standard_Monomial_Vectors_io;      use Standard_Monomial_Vectors_io;
with DoblDobl_Monomial_Vectors;
with DoblDobl_Monomial_Vectors_io;      use DoblDobl_Monomial_Vectors_io;
with QuadDobl_Monomial_Vectors;
with QuadDobl_Monomial_Vectors_io;      use QuadDobl_Monomial_Vectors_io;
with Random_Monomial_Vectors;           use Random_Monomial_Vectors;

procedure ts_monvec is

-- DESCRIPTION :
--   Tests the operations on monomial vectors.

  procedure Standard_Eval
              ( m : in Standard_Monomial_Vectors.Link_to_Monomial_Vector ) is

  -- DESCRIPTION :
  --   Evaluates m at a random vector, in double precision.

    dim : constant integer32 := integer32(m(m'first).dim);
    x : constant Standard_Complex_Vectors.Vector(1..dim)
      := Standard_Random_Vectors.Random_Vector(1,dim);
    y,z : Standard_Complex_Numbers.Complex_Number;
    yd,zd : Standard_Complex_Vectors.Vector(1..dim);
    nrm : double_float;

  begin
    y := Standard_Monomial_Vectors.Eval(m,x);
    put("y : "); put(y); new_line;
    Standard_Monomial_Vectors.Diff(m,x,yd);
    Standard_Monomial_Vectors.Speel(m,x,z,zd);
    put("z : "); put(z); new_line;
    Standard_Complex_Numbers.Sub(y,z);
    nrm := Standard_Complex_Numbers.AbsVal(y);
    put("Norm of the difference :"); put(nrm,3); new_line;
    put_line("The derivatives with the straighforward algorithm :");
    put_line(yd);
    put_line("The derivatives with the Speelpenning algorithm :");
    put_line(zd);
    Standard_Complex_Vectors.Sub(yd,zd);
    nrm := Standard_Complex_Vector_Norms.Max_Norm(yd);
    put("Max norm of the difference :"); put(nrm,3); new_line;
  end Standard_Eval;

  procedure DoblDobl_Eval
              ( m : in DoblDobl_Monomial_Vectors.Link_to_Monomial_Vector ) is

  -- DESCRIPTION :
  --   Evaluates m at a random vector, in double double precision.

    dim : constant integer32 := integer32(m(m'first).dim);
    x : constant DoblDobl_Complex_Vectors.Vector(1..dim)
      := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    y,z : DoblDobl_Complex_Numbers.Complex_Number;
    yd,zd : DoblDobl_Complex_Vectors.Vector(1..dim);
    nrm : double_float;

  begin
    y := DoblDobl_Monomial_Vectors.Eval(m,x);
    put("y : "); put(y); new_line;
    DoblDobl_Monomial_Vectors.Diff(m,x,yd);
    DoblDobl_Monomial_Vectors.Speel(m,x,z,zd);
    put("z : "); put(z); new_line;
    DoblDobl_Complex_Numbers.Sub(y,z);
    nrm := hi_part(DoblDobl_Complex_Numbers.AbsVal(y));
    put("Norm of the difference :"); put(nrm,3); new_line;
    put_line("The derivatives with the straighforward algorithm :");
    put_line(yd);
    put_line("The derivatives with the Speelpenning algorithm :");
    put_line(zd);
    DoblDobl_Complex_Vectors.Sub(yd,zd);
    nrm := hi_part(DoblDobl_Complex_Vector_Norms.Max_Norm(yd));
    put("Max norm of the difference :"); put(nrm,3); new_line;
  end DoblDobl_Eval;

  procedure QuadDobl_Eval
              ( m : in QuadDobl_Monomial_Vectors.Link_to_Monomial_Vector ) is

  -- DESCRIPTION :
  --   Evaluates m at a random vector, in quad double precision.

    dim : constant integer32 := integer32(m(m'first).dim);
    x : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    y,z : QuadDobl_Complex_Numbers.Complex_Number;
    yd,zd : QuadDobl_Complex_Vectors.Vector(1..dim);
    nrm : double_float;

  begin
    y := QuadDobl_Monomial_Vectors.Eval(m,x);
    put("y : "); put(y); new_line;
    QuadDobl_Monomial_Vectors.Diff(m,x,yd);
    QuadDobl_Monomial_Vectors.Speel(m,x,z,zd);
    put("z : "); put(z); new_line;
    QuadDobl_Complex_Numbers.Sub(y,z);
    nrm := hihi_part(QuadDobl_Complex_Numbers.AbsVal(y));
    put("Norm of the difference :"); put(nrm,3); new_line;
    put_line("The derivatives with the straighforward algorithm :");
    put_line(yd);
    put_line("The derivatives with the Speelpenning algorithm :");
    put_line(zd);
    QuadDobl_Complex_Vectors.Sub(yd,zd);
    nrm := hihi_part(QuadDobl_Complex_Vector_Norms.Max_Norm(yd));
    put("Max norm of the difference :"); put(nrm,3); new_line;
  end QuadDobl_Eval;

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Prompts for dimension, largest exponent, and then generates
  --   a random monomial vector of monomials with coefficients
  --   in standard double precision.

    size,dim : integer32 := 0;
    expmax : natural32 := 0;
    v : Standard_Monomial_Vectors.Link_to_Monomial_Vector;

  begin
    put_line("Testing monomial operations in double precision ...");
    put("Give the size of the vector : "); get(size);
    put("Give the dimension : "); get(dim);
    put("Give the largest exponent : "); get(expmax);
    declare
      ranvec : constant Standard_Monomial_Vectors.Monomial_Vector
             := Standard_Random_Monomial_Vector(size,dim,expmax,true);
    begin
      v := new Standard_Monomial_Vectors.Monomial_Vector'(ranvec);
      put_line("a random monomial vector : "); put(v);
      Standard_Eval(v);
    end;
  end Standard_Main;

  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Prompts for dimension, largest exponent, and then generates
  --   a random monomial vector of monomials with coefficients
  --   in double double precision.

    size,dim : integer32 := 0;
    expmax : natural32 := 0;
    v : DoblDobl_Monomial_Vectors.Link_to_Monomial_Vector;

  begin
    put_line("Testing monomial operations in double double precision ...");
    put("Give the size of the vector : "); get(size);
    put("Give the dimension : "); get(dim);
    put("Give the largest exponent : "); get(expmax);
    declare
      ranvec : constant DoblDobl_Monomial_Vectors.Monomial_Vector
             := DoblDobl_Random_Monomial_Vector(size,dim,expmax,true);
    begin
      v := new DoblDobl_Monomial_Vectors.Monomial_Vector'(ranvec);
      put_line("a random monomial vector : "); put(v);
      DoblDobl_Eval(v);
    end;
  end DoblDobl_Main;

  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Prompts for dimension, largest exponent, and then generates
  --   a random monomial vector of monomials with coefficients
  --   in quad double precision.

    size,dim : integer32 := 0;
    expmax : natural32 := 0;
    v : QuadDobl_Monomial_Vectors.Link_to_Monomial_Vector;

  begin
    put_line("Testing monomial operations in quad double precision ...");
    put("Give the size of the vector : "); get(size);
    put("Give the dimension : "); get(dim);
    put("Give the largest exponent : "); get(expmax);
    declare
      ranvec : constant QuadDobl_Monomial_Vectors.Monomial_Vector
             := QuadDobl_Random_Monomial_Vector(size,dim,expmax,true);
    begin
      v := new QuadDobl_Monomial_Vectors.Monomial_Vector'(ranvec);
      put_line("a random monomial vector :"); put(v);
      QuadDobl_Eval(v);
    end;
  end QuadDobl_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the precision and then lauches the test.

    ans : character;

  begin
    new_line;
    put_line("MENU to select the precision :");
    put_line("  0. double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    new_line;
    case ans is
      when '0' => Standard_Main;
      when '1' => DoblDobl_Main;
      when '2' => QuadDobl_Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_monvec;
