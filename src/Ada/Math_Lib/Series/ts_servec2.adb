with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Random_Vectors;
with Standard_Dense_Series2;
with Standard_Dense_Series2_io;          use Standard_Dense_Series2_io;
with Standard_Dense_Series2_Vectors;
with Standard_Dense_Series2_Vectors_io;  use Standard_Dense_Series2_Vectors_io;
with Standard_Dense_Series2_VecVecs;
with Standard_Dense_Vector_Series2;
with Standard_Dense_Vector_Series2_io;   use Standard_Dense_Vector_Series2_io;
with Standard_Series_Vector_Norms2;
with Random_Series_Vectors;

procedure ts_servec2 is

-- DESCRIPTION :
--   Test on vectors of truncated power series.

  procedure Write ( v : in Standard_Dense_Series2_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Writes the components of the vector to standard output.
 
  begin
    for i in v'range loop
      put("Component "); put(i,1); put_line(" :");
      put(v(i));
    end loop;
  end Write;

  procedure Standard_Test_Norm 
              ( v : in Standard_Dense_Series2_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Test on normalization of the vector v.

    use Standard_Dense_Series2;
    use Standard_Dense_Series2_Vectors;
    use Standard_Series_Vector_Norms2;

    sn : constant Series := Norm(v);
    snv : Series(sn.deg);
    nv : Vector(v'range);

  begin
    new_line;
    put_line("The norm of v : "); put(sn);
    nv := Normalize(v);
    put_line("The normalized vector : "); Write(nv);
    snv := Square_of_Norm(nv);
    put_line("Square of norm of normalized vector :"); put(snv);
  end Standard_Test_Norm;

  procedure Standard_Test ( dim,degree : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random vector of range 1..dim,
  --   with series of the given degree, in standard double precision.
  --   Tests the computation of the norm and the normalization.

    sv : Standard_Dense_Series2_Vectors.Vector(1..dim)
       := Random_Series_Vectors.Random_Series_Vector(1,dim,degree);
    vs : Standard_Dense_Vector_Series2.Vector(degree)
       := Standard_Dense_Vector_Series2.Create(sv);

  begin
    Write(sv);
    put_line("The coefficients of the vector series :"); put(vs);
    Standard_Test_Norm(sv);
  end Standard_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a dimension, and degree
  --   and then generates a random vector.

    dim,degree : integer32 := 0;

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    put("Give the degree : "); get(degree);
    Standard_Test(dim,degree);
  end Main;

begin
  Main;
end ts_servec2;
