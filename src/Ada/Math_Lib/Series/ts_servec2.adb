with text_io;                             use text_io;
with Communications_with_User;            use Communications_with_User;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;            use Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Random_Vectors;
with Standard_Dense_Series2;
with Standard_Dense_Series2_io;           use Standard_Dense_Series2_io;
with Standard_Dense_Series2_Vectors;
with Standard_Dense_Series2_Vectors_io;   use Standard_Dense_Series2_Vectors_io;
with Standard_Dense_Vector_Series2;

procedure ts_servec2 is

-- DESCRIPTION :
--   Test on vectors of truncated power series.

  function Random_Series
             ( degree : integer32 )
             return Standard_Dense_Series2.Series is

  -- DESCRIPTION :
  --   Returns a power series of the given degree 
  --   with random complex coefficients.

    cff : constant Standard_Complex_Vectors.Vector(0..degree)
        := Standard_Random_Vectors.Random_Vector(0,degree);

  begin
    return Standard_Dense_Series2.Create(cff);
  end Random_Series;

  function Random_Series_Vector
             ( first,last,degree : integer32 )
             return Standard_Dense_Series2_Vectors.Vector is

    res : Standard_Dense_Series2_Vectors.Vector(first..last);

  begin
    for k in res'range loop
      res(k) := new Standard_Dense_Series2.Series'(Random_Series(degree));
    end loop;
    return res;
  end Random_Series_Vector;

  procedure Write ( v : in Standard_Dense_Series2_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Writes the components of the vector to standard output.
 
  begin
    for i in v'range loop
      put("Component "); put(i,1); put_line(" :");
      put(v(i));
    end loop;
  end Write;

  procedure Standard_Test ( dim,degree : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random vector of range 1..dim,
  --   with series of the given degree, in standard double precision.
  --   Tests the computation of the norm and the normalization.

    sv : Standard_Dense_Series2_Vectors.Vector(1..dim)
       := Random_Series_Vector(1,dim,degree);

  begin
    Write(sv);
   -- put_line("The coefficients of the vector series :"); put(vs);
   -- Standard_Test_Norm(sv);
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
