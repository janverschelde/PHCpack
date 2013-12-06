with text_io,integer_io;                use text_io,integer_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;

procedure ts_squeeze is

-- DESCRIPTION :
--   This programs tests squeezing an n-dimensional system into the
--   plane, hoping to preserve the information about the components.

  function Random_Sum ( p : Poly_Sys ) return Poly is

  -- DESCRIPTION :
  --   Returns a random sum of the polynomials in the system.

    res : Poly := Null_Poly;
    acc : Poly;

  begin
    for i in p'range loop
      acc := Random1*p(i);
      Add(res,acc);
      Clear(acc);
    end loop;
    return res;
  end Random_Sum;

  function Random_Sums ( p : Poly_Sys ) return Poly_Sys is

  -- DESCRIPTION :
  --   Returns a polynomial system of the same range as p
  --   with random sums of the polynomials in p.

    res : Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Random_Sum(p);
    end loop;
    return res;
  end Random_Sums;

  function Product ( p : Poly_Sys ) return Poly is

  -- DESCRIPTION :
  --   Returns the product of the polynomials in p.

    res : Poly := p(p'first);

  begin
    for i in p'first+1..p'last loop
      Mul(res,p(i));
    end loop;
    return res;
  end Product;

  function Trunc2 ( p : Poly ) return Poly is

  -- DESCRIPTION :
  --   Returns a polynomial in two variables by truncation of
  --   all other variables.

    res : Poly := Null_Poly;
    d : Degrees := new Standard_Natural_Vectors.Vector(1..2);

    procedure Trunc_Term ( t : in Term; continue : out boolean ) is

      tt : Term;

    begin
      tt.dg := d;
      tt.dg(1) := t.dg(1);
      tt.dg(2) := t.dg(2);
      tt.cf := t.cf;
      Add(res,tt);
      continue := true;
    end Trunc_Term;
    procedure Trunc_Terms is new Visiting_Iterator(Trunc_Term);

  begin
    Trunc_Terms(p);
    Clear(d);
    return res;
  end Trunc2;

  procedure Squeeze ( file : in file_type; p : in Poly_Sys ) is

  -- DESCRIPTION :
  --   Makes two products of random sums of the polynomials in p and
  --   truncates the two polynomials to two variables.
  --   The squeezed system is then written to the file.

    squeezed : Poly_Sys(1..2);
    ran_sum1 : Poly_Sys(p'range) := Random_Sums(p);
    prod1 : Poly := Product(ran_sum1);
    ran_sum2 : Poly_Sys(p'range) := Random_Sums(p);
    prod2 : Poly := Product(ran_sum2);

  begin
    squeezed(1) := prod1; --Trunc2(prod1);
    squeezed(2) := prod2; --Trunc2(prod2);
    put_line(file,squeezed);
  end Squeeze;

  procedure Main is

    file : file_type;
    lp : Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Reading the polynomial system to be squeezed.");
    get(lp);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    Squeeze(file,lp.all);
  end Main;

begin
  Main;
end ts_squeeze;
